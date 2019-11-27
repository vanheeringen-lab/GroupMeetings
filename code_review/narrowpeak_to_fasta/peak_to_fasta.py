"""
From narrowpeaks -> fasta.
"""
import multiprocessing

import pybedtools
import pyfaidx


LINEWIDTH = 80
PEAKWIDTH = 100


def get_open(assembly: str, stage: str) -> str:
    """
    Convert from a narrowpeak file for the peaks of an assembly and developmental stage a fasta file
    which contains a sequence of a range of PEAKWIDTH around the summit of each peak.

    Example output fasta file:
    >lcl|chr1|peak_0
    gttaatggcgcttggcaggccgatttatatggggcattcccgccacctattggacaggagtgtgaaccgcaCGTGTTATA
    CGTTTTGCTAGCAATATGTC

    >lcl|chr2|peak_1
    TTCTCACATGATGTTGCCACTGGAAGTCGCCATGATGGTTCCCTGTGAAGAAAGCTTATCAAGACATTACTAATAGATAG
    CCGTCGAGCGTAATATTTAG

    :return: path to the fasta file
    """
    # load the chromosome
    genome = pyfaidx.Fasta(f"/vol/genomes/{assembly}/{assembly}.fa")

    # load the peak file
    narrowpeaks = pybedtools.BedTool(f"/vol/atac-seq/genrich/{assembly}-{stage}_peaks.narrowPeak")

    fasta_strings = []
    for narrowpeak in narrowpeaks:
        if narrowpeak.length < 2000:
            # get the summit and the flanking low and high sequences
            summit = narrowpeak.start + int(narrowpeak.fields[-1])
            low, high = summit - PEAKWIDTH // 2, summit + PEAKWIDTH // 2

            # ignore if out of bounds
            if low < 0 or high > len(genome[narrowpeak.chrom]):
                continue

            # get the sequence and split into chunks of max line width
            sequence = str(genome[narrowpeak.chrom][low:high])
            sequence = '\n'.join(sequence[i:i + LINEWIDTH]
                                 for i in range(0, len(sequence), LINEWIDTH))

            # store the sequence as fasta (NCBI local)
            fasta_strings.append(f">lcl|{narrowpeak.chrom}|{narrowpeak.name}\n"
                                 f"{sequence}")

    # combine all the sequences into a single fasta and save
    fasta = '\n\n'.join(fasta_strings)
    with open(f"/vol/peaks/{assembly}-{stage}_open.fa", "w") as output_fasta:
        output_fasta.write(fasta)
    return f"/vol/peaks/{assembly}-{stage}_open.fa"


def get_closed(assembly: str, stages: list) -> None:
    """
    Convert from a list of narrowpeak files for the peaks of an assembly across developmental stages
    a fasta file which contains a sequence of a range of PEAKWIDTH around randomly chosen areas not
    in any of the narrowpeak files.

    Example output fasta file:
    >lcl|chr10|9650
    CCAAGCAGCAGCACTCGGGCACATTCATTATGTAAGGACGACAAGCCCCGCCCACATTAACACCCCCCCCCCCCATCCTA
    GTTCCACCCACATTGACACC

    >lcl|chr15|32389
    ACAAACAAAACATTAATGGAATGAGTATAATAATTGTAACATTCTTAATCGATCATAACTTCCTTTAAGAGGAAGACGAT
    ATGTTTTTGTCTTAACGTAC

    :return: path to the fasta file
    """
    # load the chromosome
    genome = pyfaidx.Fasta(f"/vol/genomes/{assembly}/{assembly}.fa")

    # load the peak files
    nr_peaks = 0
    combination = pybedtools.BedTool(())
    for stage in stages:
        bed = pybedtools.BedTool(f"/vol/atac-seq/genrich/{assembly}-{stage}_peaks.narrowPeak")
        nr_peaks += len(bed)
        combination = combination.cat(bed)

    combination = combination.slop(
        l=PEAKWIDTH,
        r=0,
        g=f"/vol/genomes/{assembly}/{assembly}.fa.sizes"
    ).merge()

    # store
    closed_sizes = []
    closed_fasta = []
    for i, closed_region in enumerate(combination):
        sequence = str(genome[closed_region.chrom][closed_region.start:closed_region.end])
        sequence = '\n'.join(sequence[i:i + LINEWIDTH]
                             for i in range(0, len(sequence), LINEWIDTH))

        length = closed_region.end - closed_region.start
        closed_sizes.append(f"{closed_region.chrom}_{i}\t{length}\t{closed_region.start}")
        closed_fasta.append(f">{closed_region.chrom}_{i}\n"
                            f"{sequence}")

    #
    with open(f"/vol/peaks/{assembly}_closed.fa.sizes", 'w') as sizes:
        sizes.write('\n'.join(closed_sizes))
    with open(f"/vol/peaks/{assembly}_closed.fa", 'w') as fasta:
        fasta.write('\n'.join(closed_fasta))
    genome = pyfaidx.Fasta(f"/vol/peaks/{assembly}_closed.fa")

    poss = combination.random(l=PEAKWIDTH,
                              n=nr_peaks,
                              g=f"/vol/peaks/{assembly}_closed.fa.sizes")

    fasta_strings = []
    for pos in poss:
        sequence = str(genome[pos.chrom][pos.start:pos.end])
        sequence = '\n'.join(sequence[i:i + LINEWIDTH]
                             for i in range(0, len(sequence), LINEWIDTH))

        # store the sequence as fasta (NCBI local)
        chrom = '_'.join(pos.chrom.split('_')[:-1])
        fasta_strings.append(f">lcl|{chrom}|{pos.chrom.split('_')[-1]}\n"
                             f"{sequence}")

    fasta = '\n\n'.join(fasta_strings)
    with open(f"/vol/peaks/{assembly}_closed.fa", "w") as output_fasta:
        output_fasta.write(fasta)


data = {'danRer11': ['48h', '8somites', '80%epiboly', 'dome', 'shield'],
        'mm10':     ['early_2cell', '2cell', '4cell', '8cell',
                     'Mm_E7', 'Mm_E8', 'Mm_E9', 'Mm_E10', 'Mm_E12', 'Mm_E14', 'Mm_E16', 'Mm_E18'],
        'Spur_5.0': ['28h', '18h', '24h', '30h', '39h', '50h', '60h', '70h'],
        'oryLat2':  ['st11', 'st13', 'st15', 'st19', 'st21', 'st24', 'st25', 'st28',
                     'st32', 'st36', 'st40']}

pool = multiprocessing.Pool(25)
pool.starmap(get_closed, (data.items()))
pool.starmap(get_open, [(assembly, stage)
                        for assembly, stages in data.items() for stage in stages])
