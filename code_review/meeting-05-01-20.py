from app import db, login


class Base_Class():
    """ Base class, has a usefull function to set data.
        Is inherited by all other database objects. """

    def set_data(self, **kwargs):
        """ Tries to find all attributes that correspong to the given arguments.
            If the attributes exist the data is set. Accepts a dictionary. """

        for attr, value in kwargs.items():
            # Attribute is found and data is set.
            if hasattr(self, attr):
                setattr(self, attr, value)

            # No attribute found.
            else:
                print(f'{self} has no attribute called "{attr}"')


class Flowcell(Base_Class, db.Model):
    """ Flowcell class. Besides standard database data,
        it contains functions to detect barcode clashes. """

    __tablename__ = 'flowcell'
    ID = db.Column(db.Integer, primary_key=True)

    #  Foreign keys and relations
    ID_operator = db.Column(db.Integer, db.ForeignKey('user.ID'))
    operator = db.relationship('User', backref='flowcells',
                               foreign_keys=[ID_operator], lazy=True)

    #  Data
    name = db.Column(db.String(128), index=True)
    buffer = db.Column(db.Integer)
    reagents = db.Column(db.Integer)
    run_type = db.Column(db.String(64))
    barcode_length = db.Column(db.Integer)

    def __repr__(self):
        return f'{self.name} ({self.ID})'

    def compare_samples(self, samples):
        """ Compare all barcodes attached to samples, accepts a list. """
        clashes = []

        for n, sample_1 in enumerate(samples[:-1]):
            for sample_2 in samples[n+1:]:
                # Sample 1 vs sample 2
                # Each sample has a variable called barcodes,
                # which is a list containing multiple Barcode objects
                for barcode_1 in sample_1.barcodes:
                    for barcode_2 in sample_2.barcodes:
                        # Compare sequence(s) belonging to the barcodes.
                        if not compare_barcodes(barcode_1.sequences,
                                                barcode_2.sequences):
                            # Append clashing samples and barcodes.
                            # String is this way because of limited PEP8
                            # line length....
                            clashes.append(f'{sample_1} ({barcode_1}) vs' +
                                           f'{sample_2} ({barcode2})')

        # Return the clashes to routes.py, if it's empty there's no problem.
        # Otherwise the clashes can be displayed.
        return clashes


def compare_barcodes(first, second):
    """ Compare two barcodes, accepts two strings. """

    # Sort by length, shortest sequence is needed for ratio calculation.
    shortest, longest = sorted([first, second], key=len)

    # Map the two sequences. Map as a 1 if exact same nucleotide,
    # otherwise map a 0.
    result = map(lambda x, y: 1 if x == y else 0,
                 shortest, longest)

    # Calculate percentage of similarity.
    similarity = result.count(1) / len(shortest)

    if similarity > 0.76:
        # 0.76 seems to give the same clash results as Rita's list.
        return False
    else:
        # No clash.
        return True
