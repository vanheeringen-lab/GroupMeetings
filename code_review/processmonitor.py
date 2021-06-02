#!/usr/bin/env python3
import os
import grp
import psutil
import subprocess
import datetime
from socket import gethostname

MAIL_WHITELIST_GROUPS = ["slrinzema"]
RAM_CUTOFF = 0.01
RAM_TOTAL_CUTOFF = 10.0
HOST = gethostname()
LOG_ROOT = f"/ceph/rimlsfnwi/data/moldevbio/heeringen"\
    f"/slrinzema/sleeplogs/{HOST}"
LOG_FILE = ""


def main():
    """ Loop over current processes and seperate them by user.
        Afterwards write log file & send emails if necessary
    """

    _init_log()

    # Sort processes by users
    users = {}  # keys: username, values: User object.
    for p in psutil.process_iter():
        u = p.username()
        if u == "root":
            continue
        if u not in users:
            users[u] = User(u)
        users[u].add_process(p)

    # Write the logs and send emails
    with open(LOG_FILE, "w") as log:
        # Sort users by total memory used
        sorted_users = sorted(users.values(),
                              key=lambda x: x.total_memory,
                              reverse=True)

        for user in sorted_users:
            # Send mail if needed
            _send_email(user)

            # write into log
            log.write("="*80 + "\n")
            log.write(str(user))


def _init_log():
    # Initialize a log folder by date & a log file name
    global LOG_ROOT, LOG_FILE

    now = datetime.datetime.now()
    today = now.strftime("%Y-%m-%d")
    log_dir = os.path.join(LOG_ROOT, today)
    # Make necessary directories
    if not os.path.isdir(LOG_ROOT):
        os.mkdir(LOG_ROOT)
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)

    # Make log file name by current time
    LOG_FILE = os.path.join(log_dir, now.strftime("%H:%M.log"))


def _user_whitelisted(user):
    # Check if user is on whitelist
    if user in MAIL_WHITELIST_GROUPS:
        return True

    groups = [g.gr_name for g in grp.getgrall() if user in g.gr_mem]
    return any(g in groups for g in MAIL_WHITELIST_GROUPS)


def _send_email(user):
    # Send email if user is on whitelist
    if not _user_whitelisted(user) \
            or user.total_memory_sleeping < RAM_TOTAL_CUTOFF:
        return

    message = f"Hello {user.name},\n\n" + \
              f"On {HOST} you currently have {len(user.processes_sleeping)}" + \
              f"processes that are sleeping and using" + \
              f"{user.total_memory_sleeping}% of total memory.\n" + \
              f"More info:\n\n{user}"
    command = f"echo '{message}' |" + \
              f"mail -s 'Sleeping processes on {HOST}' {user.name}"
    subprocess.Popen(command, shell=True)


class User:
    """ A class to hold all process data for a user. """

    def __init__(self, name):
        """ Init user with name. """
        self.name = name
        self.processes = {}


    def add_process(self, process):
        """ Add a process to the processes dictionary. """
        self.processes[process.pid] = process


    @property
    def processes_total(self):
        """ Get total processes. """
        return len(self.processes)


    @property
    def processes_by_cutoff(self):
        """ Seperate processes by ram usage.

            Returns two lists 'above' and 'below'.
        """
        above = []
        below = []
        for pid, process in self.processes.items():
            if psutil.pid_exists(pid) \
                    and process.memory_percent() <= RAM_CUTOFF:
                below.append(process)
            else:
                above.append(process)
        return above, below


    @property
    def processes_sleeping(self):
        """ Gets all sleeping processes and returns them as a list. """
        sleeping = []
        for pid, process in self.processes.items():
            if psutil.pid_exists(pid) and \
                    process.status() == "sleeping":
                sleeping.append(process)
        return sleeping


    @property
    def total_memory(self):
        """ Get memory usage of all processes. """
        mem = 0
        for pid, process in self.processes.items():
            if psutil.pid_exists(pid):
                mem += process.memory_percent()
        return round(mem, 3)


    @property
    def total_memory_sleeping(self):
        """ Get memory usage of all sleeping processes. """
        mem = 0
        for pid, process in self.processes.items():
            if psutil.pid_exists(pid) and process.status() == "sleeping":
                mem += process.memory_percent()
        return round(mem, 3)


    def _process_info(self, process):
        # Get info for the processes and turn it into a readable string
        startdate = datetime.datetime.fromtimestamp(
                process.create_time()).strftime("%Y-%m-%d %H:%M:%S")
        repr = f"{process.pid}\t{process.name()}\t{round(process.memory_percent(),2)}%\t{startdate}\t{process.status()}\n"
        cmdline = ' '.join(process.cmdline())
        repr += f" └─ {cmdline}\n"
        return repr


    def __repr__(self):
        # Return printable info for this object.
        above, below = self.processes_by_cutoff
        repr = f"{self.name}\n{self.processes_total} processes, {len(above)} above {RAM_CUTOFF}% RAM\n{self.total_memory}% RAM total\n"
        if self.processes_total == 0:
            return repr

        if len(above) > 0:
            repr += f"\nProcesses above {RAM_CUTOFF}% RAM:\n"
            repr += "pid\tname\tmem %\tstartdate\tstate\n"
            for process in above:
                repr += self._process_info(process)

        if len(below) > 0:
            repr += f"\nProcesses below {RAM_CUTOFF}% RAM:\n"
            repr += "pid\tname\tmem %\tstartdate\tstate\n"
            for process in below:
                repr += self._process_info(process)

        return repr + "\n"


if __name__ == "__main__":
    main()
