# Progress/struck

"""
MultiQC gets stuck if it needs to parse a large amount of samples
Bash also crashes after ~100 arguments.

Ideally, rule MultiQC can see the need to split its input in multiple batches,
Input files then need to be split by sample name, but these can change
between steps by replicate merging.

Below it the current QC file aggregation code
from https://github.com/vanheeringen-lab/snakemake-workflows/blob/replicates/rules/qc.smk
"""

def get_qc_files(wildcards):
    assert 'quality_control' in globals(), "When trying to generate multiqc output, make sure that the "\
                                           "variable 'quality_control' exists and contains all the "\
                                           "relevant quality control functions."
    qc = []
    if 'replicate' in samples and config.get('technical_replicates') == 'merge':
        # trimming qc on individual samples
        for sample in samples[samples['assembly'] == wildcards.assembly].index:
            qc.extend(get_trimming_qc(sample))

        # qc on merged technical replicates
        for replicate in set(samples[samples['assembly'] == wildcards.assembly].replicate):
            for function in [func for func in quality_control if
                             func.__name__ not in ['get_peak_calling_qc', 'get_trimming_qc']]:
                qc.extend(function(replicate))

        # qc on combined biological replicates
        if get_peak_calling_qc in qc:
            for condition in set(samples[samples['assembly'] == wildcards.assembly].condition):
                qc.extend(function(condition))
    else:
        for sample in samples[samples['assembly'] == wildcards.assembly].index:
            for function in quality_control:
                qc.extend(function(sample))
    return qc




# Code review-ish part:

"""
Lets install a nice and small genome:
    genomepy install ASM2732v1 NCBI -a


ASM2732v1.fa has only one contig:

>ANONYMOUS NC_000908.2 Mycoplasma genitalium G37, complete sequence
TAAGTTATTATTTAGTTAATACTTTTAACAATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATAGTTAA
ATACCTTCCTTAATACTGTTAAATTATATTCAATCAATACATATATAATATTATTAAAATACTTGATAAGTATTATTTAG
ATATTAGACAAATACTAATTTTATATTGCTTTAATACTTAATAAATACTACTTATGTATTAAGTAAATATTACTGTAATA


ASM2732v1.annotation.gtf.gz looks a little different:

NC_000908.2	/home/siebrenf/.local/share/genomes/ASM2732v1/tmpkk6wfo7o/ASM2732v1.annotation.gp	transcript	579224	580033	.	-	.	gene_id "gene-MG_RS02770"; transcript_id "cds-WP_009885563.1";  gene_name "gene-MG_RS02770";
NC_000908.2	/home/siebrenf/.local/share/genomes/ASM2732v1/tmpkk6wfo7o/ASM2732v1.annotation.gp	exon	579224	580033	.	-	.	gene_id "gene-MG_RS02770"; transcript_id "cds-WP_009885563.1"; exon_number "1"; exon_id "cds-WP_009885563.1.1"; gene_name "gene-MG_RS02770";
NC_000908.2	/home/siebrenf/.local/share/genomes/ASM2732v1/tmpkk6wfo7o/ASM2732v1.annotation.gp	CDS	579227	580033	.	-	0	gene_id "gene-MG_RS02770"; transcript_id "cds-WP_009885563.1"; exon_number "1"; exon_id "cds-WP_009885563.1.1"; gene_name "gene-MG_RS02770";


Note here that the contig name in the GTF matches the SECOND name in the fasta.
The official fasta format dictates the FIRST name should match between fasta and GTF.
Tools like STAR use this format, and therefore crash.

Since *most* fastas have the matching contig name (somewhere) in their header as well, we can try to fix this.
That is what the following function does...

"""

import os.path
import sys
import subprocess as sp

from tempfile import TemporaryDirectory
from genomepy.utils import bgunzip_and_name, bgrezip


def sanitize_annotation(genome, gtf_file=None, sizes_file=None, out_dir=None):
    """
    Matches the toplevel sequence names in annotation.gtf to those in genome.fa.

    The fasta and gtf formats dictate the 1st field in the genome header and gtf should match.
    In some cases the genome.fa has multiple names per header, the 1st field not matching those in the gtf.
    If this occurs a conversion table can be made to rename the sequence names in either file.

    This script changes the names in the gtf to match those in the genome.fa, if needed.

    genome: class object
        Genome class object (contains genome name and file paths)

    sizes_file: str , optional
        Path to the genome sizes file

    gtf_file: str , optional
        Path to the (gzipped) annotation.gtf file

    out_dir: str , optional
        Location of the genome (and the README.txt)
    """
    # set default file paths if none are given
    if out_dir is None:
        out_dir = os.path.dirname(genome.filename)
    if gtf_file is None:
        gtf_file = os.path.join(out_dir, genome.name + ".annotation.gtf.gz")
    if sizes_file is None:
        sizes_file = os.path.join(out_dir, genome.name + ".fa.sizes")

    # unzip gtf if zipped and return up-to-date name
    if gtf_file.endswith(".gz"):
        sp.check_call("gunzip -f {}".format(gtf_file), shell=True)
        gtf_file = gtf_file[:-3]

    # Check if genome and annotation have matching chromosome/scaffold names
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                gtf_id = line.split("\t")[0]
                break

    with open(sizes_file, "r") as sizes:
        for line in sizes:
            fa_id = line.split("\t")[0]
            if fa_id == gtf_id:
                sp.check_call("gzip -f {}".format(gtf_file), shell=True)

                readme = os.path.join(out_dir, "README.txt")
                with open(readme, "a") as f:
                    f.write("genome and annotation have matching sequence names.\n")
                return

    # generate a gtf with matching scaffold/chromosome IDs
    sys.stderr.write(
        "Genome and annotation do not have matching sequence names! Creating matching annotation files...\n"
    )

    # unzip genome if zipped and return up-to-date genome name
    bgzip, genome_file = bgunzip_and_name(genome)

    # determine which element in the fasta header contains the location identifiers used in the annotation.gtf
    header = []
    with open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                header = line.strip(">\n").split(" ")
                break

    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                loc_id = line.strip().split("\t")[0]
                try:
                    element = header.index(loc_id)
                    break
                except ValueError:
                    continue
        else:
            # re-zip genome if unzipped
            bgrezip(bgzip, genome_file)

            sys.stderr.write(
                "WARNING: Cannot correct annotation files automatically!\n"
                + "Leaving original version in place."
            )
            return

    # build a conversion table
    ids = {}
    with open(genome_file, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                line = line.strip(">\n").split(" ")
                if line[element] not in ids.keys():
                    ids.update({line[element]: line[0]})

    # re-zip genome if unzipped
    bgrezip(bgzip, genome_file)

    # remove old files
    with TemporaryDirectory(dir=out_dir) as tmpdir:
        old_gtf_file = os.path.join(tmpdir, os.path.basename(gtf_file))
        bed_file = gtf_file.replace("gtf", "bed")
        sp.check_call("mv {} {}".format(gtf_file, old_gtf_file), shell=True)
        sp.check_call("rm {}".format(bed_file + ".gz"), shell=True)

        # rename the location identifier in the new gtf (using the conversion table)
        with open(old_gtf_file, "r") as oldgtf, open(gtf_file, "w") as newgtf:
            for line in oldgtf:
                line = line.split("\t")
                line[0] = ids[line[0]]
                line = "\t".join(line)
                newgtf.write(line)

    # generate a bed with matching scaffold/chromosome IDs
    cmd = (
        "gtfToGenePred {0} /dev/stdout | " "genePredToBed /dev/stdin {1} && gzip -f {1}"
    )
    sp.check_call(cmd.format(gtf_file, bed_file), shell=True)
    sp.check_call("gzip -f {}".format(gtf_file), shell=True)

    readme = os.path.join(out_dir, "README.txt")
    with open(readme, "a") as f:
        f.write("corrected annotation files generated succesfully.\n")
    return
