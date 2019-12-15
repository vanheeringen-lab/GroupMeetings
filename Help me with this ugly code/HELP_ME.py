@app.route('/view/submitted',   methods=['GET', 'POST'])
@login_required
def viewSubmitted():
    """ Function to prepare the site view containing all submitted samples.
        It checks for barcode clashes and makes a pipette scheme.

        TODO: make nicer, especially the pipette scheme generating part.
              But how? """

    # Sort all samples with a state of 2 (AKA submitted)
    # by experiment name, results in nice seperated batches.
    samples = sorted(Sample.query.filter_by(state=2).all(),
                     key=operator.attrgetter('experiment.name'))

    # Make a dictionairy and store samples by experiment.
    samples_sorted = {}
    for sample in samples:
        if sample.experiment not in samples_sorted:
            samples_sorted[sample.experiment] = [sample]
        else:
            samples_sorted[sample.experiment].append(sample)

    # If the form on the site is submitted this will result in True
    if request.method == "POST":
        # Submit button has been clicked.
        # Parse form request into sample id's.
        samples = [int(x) for x in request.form.getlist('checkbox_sample')]

        # Check for barcode clashes, done with tools/compare_barcodess.py
        # Pass on id's only, had trouble passing on the "Sample" classes
        clashes = loop_samples(samples, use_ids=True)
        if clashes != {}:
            # clashes are found, make a nice error message.
            message = 'Barcode clashes\n'
            for sampleComb, clash in clashes.items():
                message += '{} {}\n'.format(sampleComb, clash)
            flash(message)  # show error message on site

        # No clashes, and also there are samples.
        elif samples != []:
            # Get Sample objects and sort by researcher
            samples = sorted([Sample.query.filter_by(id=int(x)).first()
                             for x in samples],
                             key=operator.attrgetter('researcher.fullname'))

            # TODO: replace 'testfile' with a name reflecting date
            with open('testfile.xlsx', 'w') as reqfile:
                researchers = []
                for sample in samples:
                    if researchers == []:
                        researchers.append(sample.researcher)
                        # TODO: make this nice
                        reqfile.write("""# Headed, Escapable delimited Tabular text data. Delimiter=  /(\t)/
#
#
#
#  c1:\tresearcher::    Full, normalized name of the researcher.  domainFrom=   researcher.  Format=    personName.
#  c2:\tlibNr:: library number.  inheritFrom=   sampleLibs.  uniqueWith=    -.  Format= integer.
#  c3:\tsampleName::    Sample name.  Context=  libNr.  Format= splitName.
#  c4:\tnDNAc:: DNA concentration, numeric part.
#  c5:\tsDNAc:: 0
#  c6:\tfragBp::    Central tendency of fragment size in number of base pairs.  Format= integer.  Length=   3
#  c7:\tadapterKit::    Reference to the adapter kit.  Use '-' if you want to supply barcodes instead of indices.  domainFrom=  kitBarcodes.  Format=   text.  Allowed= -.
#  c8:\tadaNrs
#  c9:\treqMR:: Requested milion reads
# c10:\tmethodType::    The cells were prepared for sequencing using method \methodType.  domainFrom=   methodType.  Allowed=   ;.
# c11:\tadaID:: Unique ID of adapter (not barcode) within department, with unique letter for the adapterKit and the number within the adapterKit.  Barcodes of adapters from different adapterKits may be not be different enough, even equal!
# c12:\tuseBarClash::   Sorted list of adaIDs that can in no way be combined with this one.  Within current standard adapterKits, only other current standard adapterKits are listed, while for old and special adapterKits, all clashes are listed.
# c13:\tseqSpec::   Specify how to sequence if that may be unclear or as a reminder (about previous or for future discussion) if it is special.  Use '-' for none.  Format= text.  Allowed= ;-.
# c14:\tDNAconc::   DNA concentration.  This is including the unit (e.g. 3 ng/µl; 3 ng/ul; 0.4 nM).  Format=    physicalValue.
#
#researcher\tlibNr\tsampleName\tnDNAc\tsDNAc\tfragBp\tadapterKit\tadaNrs\treqMR\tmethodType\tadaID\tuseBarClash\tseqSpec\tDNAconc\n""")

                    if sample.researcher not in researchers:
                        reqfile.write('\n\n\n')
                    reqfile.write(sample.researcher.fullname +
                                  '\t' + str(sample.id) +
                                  '\t' + str(sample.name) +
                                  '\t' + str(sample.DNAconc) +
                                  '\t' + str('-') +
                                  '\t' + str(sample.fragBP) +
                                  '\t' + str(sample.barcodes[0].adapterKit.name) +
                                  '\t' + str(''.join(x.barcode+';'
                                             for x in sample.barcodes)[:-1]) +
                                  '\t' + str(sample.sugMR) +
                                  '\t' + str(sample.experiment.methodType) +
                                  '\t' + str(''.join(x.id+';' for x in sample.barcodes)[:-1]) +
                                  '\t' + str('-') +
                                  '\t' + str(sample.seqSpec) +
                                  '\t' + str(sample.DNAconc) +
                                  '\n')

        return send_file('/home/slrinzema/SeqAdmin/testfile.xlsx', as_attachment=True)

    return render_template('viewSubmitted.html', title='View Submitted Samples',
                           samples=samples_sorted)


def loop_samples(samples, use_ids=False):
    """ Accepts a list of int's. Convert these to Sample objects and loop over them,
        checking for barcode clashes. Barcode clash occurs with a >75% similarity. """

    if use_ids:
        samples = getSamples(samples)  # Get Sample objects with id.

    clashes = {}  # Dictionary to hold clashes
    for N, firstSample in enumerate(samples[:-1]):
        for secondSample in samples[(N+1):]:
            # Each sample can have multiple barcodes,
            # also check out this walrus operator!
            if (found_clashes := loop_barcodes(firstSample
                                               secondSample)) != []:
                # Make a string to show on the site,
                # so people know what samples clash.
                clash_combi = "{} vs {}".format(firstSample.name,
                                                secondSample.name)
                clashes[clash_combi] = found_clashes
    return clashes


def loop_barcodes(firstSample, secondSample):
    """ Loop over all barcodes, pertaining to two samples. """

    clashes = []
    # First loop for all barcodes in first sample
    for fBarcode in firstSample.barcodes:
        # Second loop for all barcodes in second sample
        for sBarcode in secondSample.barcodes:
            if not compare_barcodes(fBarcode.barcode,
                                    sBarcode.barcode):
                # Function returned False, aka not compatible.
                barcodeClash = "{}/{}".format(fBarcode.id,
                                              sBarcode.id)
                clashes.append(barcodeClash)
    return clashes


def compare_barcodes(first, second):
    """ Compare two barcodes, returning False if they aren't compatible
        and True if they are. """

    # Sort by length, compare shortest to longest
    shortest, longest = sorted([first, second], key=len)
    similarities = 0
    for N, nuc in enumerate(shortest):
        similarities += 1 if nuc == longest[N] else 0
    ratio = similarities / len(shortest)

    if ratio > 0.76:
        # 0.76 seems to give the same results as Rita's list
        # Also if higher, no barcode compatibility.
        return False
    else:
        # Barcodes compatible
        return True
