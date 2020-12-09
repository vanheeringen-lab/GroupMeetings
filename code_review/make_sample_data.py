def make_sample_data(sample, barcode_length, barcode_amount):
    """ 
    Return a list of SampleSheet formatted data, derived from the sample() object. 
    The sample object contains info about the barcodes used. Some barcodes contain 
    sequences for both ends of the sample (i5 and i7) seperated by a + and some 
    barcodes contain multiple sequences seperated by a ;.
        
    Parameters
    ----------
    sample: object
        A Python object generated by sqlalchemy.
    barcode_length: int
        The length of barcode read by the Illumina machine.
    barcode_amount: int
        The amount of barcodes, 1 means only the i7. 2 means i7 + i5.
    
    Returns a list of strings containing samplesheet data like so:
    SampleID,SampleName,i7_index,i5_index
    """
    
    # The sample() object contains a variable called end, indicating which end
    # of the sample it attaches to. 0 means both i7 (3') and i5 (5'), 1 means i7
    # and 2 means i5. The "0" barcodes sequences are split by a + sign.
    both_sides = list(filter(lambda x: x.end == 0, sample.barcodes))

    i7_barcodes = list(filter(lambda x: x.end == 1, sample.barcodes))
    i5_barcodes = list(filter(lambda x: x.end == 2, sample.barcodes))
    for barcode in both_sides:
        i7_barcodes.append(barcode)
        i5_barcodes.append(barcode)

    # Make all possible combinations between i7 and i5 sequences.
    combinations = []
    for i7 in i7_barcodes:
        # i5 barcode not needed
        if barcode_amount == 1:
            combinations.append((i7, ""))
            continue

        # i5 barcode is needed, but the sample has none. Use the primer remnant
        # belonging to the sample's adapter kit.
        if barcode_amount == 2 and not i5_barcodes and i7.adapter_kit:
            remnant = i7.adapter_kit.primer_remnant
            combinations.append((i7, remnant))
            continue

        # Make combinations if both barcodes needed.
        for i5 in i5_barcodes:
            combinations.append((i7, i5))

    n = 0
    data = []
    for i7, i5 in combinations:
        # Sometimes a barcode objects contains multiple sequences, split by ;
        # These are always at the same end of the sample, also replaces spaces
        # just to be safe.
        i7_seq_list = i7.sequences.replace(" ", "").split(";")
        for i7_seq in i7_seq_list:
            i7_seq = adjust_sequence(i7_seq, barcode_length)
            sample_id = f"{sample.id}_{n}"
            sample_name = (f"{sample.name}"
                           # Add N for uniqueness if needed.
                           + (f"_{n}" if len(i7_seq_list) > 1 else "") 
                           + f"_{i7.name}")

            i5_seq = ""
            if barcode_amount > 1:
                if isinstance(i5, str): # This means it's a primer remnant
                    i5_seq = adjust_sequence(i5, barcode_length)
                else:
                    i5_seq = adjust_sequence(i5.sequences,
                                             barcode_length,
                                             True)
                    sample_name += f"_{i5.name}"

            sample_data = f"{sample_id},{sample_name},{i7_seq},{i5_seq}"
            data.append(sample_data)
            n += 1

    return data


def adjust_sequence(sequence, length, recomplement=False):
    """
    Supplements a sequence with AT nucleotides and returns it with the correct length.
    If recomplement is True, the sequence get reverse complemented first.
    """

    if recomplement:
        sequence = reverse_complement(sequence)
    
    # If a barcode is shorter than the read length, the machine reads AT's
    # at the empty positions.
    sequence += "AT"*20
    return sequence[:length]


def reverse_complement(seq):
    """ Return a reverse complemented sequence. """
    comp = {"A": "t", "T": "a",
            "G": "c", "C": "g"}

    seq = seq.upper()
    for key, value in comp.items():
        seq = seq.replace(key, value)

    return seq[::-1].upper()