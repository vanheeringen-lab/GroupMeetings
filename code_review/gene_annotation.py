import os
import sys
from typing import Tuple, Iterable

from loguru import logger
import mygene
from genomepy.provider import ProviderBase, cached
from genomepy import Genome
import pandas as pd


logger.remove()
logger.add(sys.stderr, format="<green>{time:YYYY-MM-DD at HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}", level="INFO")


#@cached
def ensembl_genome_info(genome_name: str) -> Tuple[str, str, str]:
    """Return Ensembl genome information for a local genome managed by genomepy.

    Parameters
    ----------
    genome_name : str
        Name of local genome.

    Returns
    -------
    (str, str, str)
        Ensembl name, accession, taxonomy_id
    """
    # Fast lookup for some common queries
    common_names = {
        "danRer11": "GRCz11",
        "hg38": "GRCh38",
        "mm10": "GRCm38",
        "dm6": "BDGP6.28",
    }
    if genome_name in common_names:
        search_term = common_names[genome_name]
    else:
        try:
            genome = Genome(genome_name)
            search_term = genome.tax_id
        except FileNotFoundError:
            logger.info(f"Genome {genome_name} not installed locally")
            p = ProviderBase.create("Ensembl")
            for name, *_rest in p.search(genome_name):
                if name == genome_name:
                    logger.info(
                        f"It can be downloaded from Ensembl: genomepy install {name} Ensembl --annotation"
                    )
                    return None
            return None

    # search Ensembl by taxonomy_id or by specific Ensembl name (if we know it)
    p = ProviderBase.create("Ensembl")
    name, accession, species, tax_id, *rest = [row for row in p.search(search_term)][0]

    # Check if the assembly_id of the current Ensembl genome is the same as the
    # local genome. If it is identical, we can correctly assume that the genomes
    # sequences are identical.
    # For the genomes in the lookup table, we already know they match.
    if genome_name in common_names or accession == genome.assembly_accession:
        return name, accession, tax_id
    else:
        print(f"Could not find a matching genome in Ensembl")
        return None


#@cached
def ncbi_assembly_report(asm_acc: str) -> pd.DataFrame:
    """Retrieve the NCBI assembly report as a DataFrame.

    Parameters
    ----------
    asm_acc : str
        Assembly accession (GCA or GCF)

    Returns
    -------
    pandas.DataFrame
        NCBI assembly report.
    """
    p = ProviderBase.create("NCBI")
    ncbi_search = list(p.search(asm_acc))
    if len(ncbi_search) > 1:
        raise Exception("More than one genome for accession")
    else:
        ncbi_name = ncbi_search[0][0].replace(" ", "_")

    # NCBI FTP location of assembly report
    logger.info(f"Found NCBI assembly {asm_acc} with name {ncbi_name}")
    assembly_report = (
        f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{asm_acc[0:3]}/"
        + f"{asm_acc[4:7]}/{asm_acc[7:10]}/{asm_acc[10:13]}/"
        + f"{asm_acc}_{ncbi_name}/{asm_acc}_{ncbi_name}_assembly_report.txt"
    )

    logger.info(f"Downloading {assembly_report}")
    header = [
        "Sequence-Name",
        "Sequence-Role",
        "Assigned-Molecule",
        "Assigned-Molecule-Location/Type",
        "GenBank-Accn",
        "Relationship",
        "RefSeq-Accn",
        "Assembly-Unit",
        "Sequence-Length",
        "UCSC-style-name",
    ]
    asm_report = pd.read_csv(assembly_report, sep="\t", comment="#", names=header)
    return asm_report


#@cached
def load_mapping(genome_name):
    logger.info("Loading chromosome mapping.")
    genome = Genome(genome_name)
    asm_acc = genome.assembly_accession

    if genome.provider not in ["UCSC", "NCBI"]:
        logger.error(f"Can't map to provider {genome.provider}")
        return None

    asm_report = ncbi_assembly_report(asm_acc)
    asm_report.loc[
        asm_report["Sequence-Role"] != "assembled-molecule", "Assigned-Molecule"
    ] = "na"

    mapping = asm_report[["Sequence-Name", "UCSC-style-name", "Assigned-Molecule", "GenBank-Accn"]]                                                                

    if genome.provider == "NCBI":
        logger.info("Mapping to NCBI sequence names")
        id_column = "Sequence-Name"
    elif genome.provider == "UCSC":
        logger.info("Mapping to UCSC sequence names")
        id_column = "UCSC-style-name"
    mapping = pd.melt(mapping, id_vars=[id_column])
    mapping = mapping[mapping["value"] != "na"]                                                                                                                
    mapping = mapping.drop_duplicates().set_index("value")[[id_column]]
    mapping.columns = ["chrom"]
    return mapping


def _local_gene_annotation(genes: Iterable[str], genome: str) -> pd.DataFrame:
    """Retrieve gene location from local annotation.

    Parameters
    ----------
    genes : Iterable
        List of gene names or gene identifiers such as ensembl_id.
    genome : str
        Genome name

    Returns
    -------
    pandas.DataFrame with gene annotation.
    """
    g = Genome(genome)
    gene_list = list(genes)
    bed = os.path.join(os.path.dirname(g.filename), f"{genome}.annotation.bed.gz")
    gene_info = pd.DataFrame()
    if os.path.exists(bed):
        df = pd.read_table(
            bed,
            index_col=3,
            usecols=[0, 1, 2, 3, 5],
            names=["chrom", "start", "end", "name", "strand"],
        )
        gene_info = df.loc[gene_list]

    # If we find more than half of the genes we assume this worked.
    if gene_info.shape[0] >= 0.5 * len(gene_list):
        return gene_info.reset_index()[["chrom", "start", "end", "name", "strand"]]


def gene_annotation(genes: Iterable[str], genome: str) -> pd.DataFrame:
    """Retrieve genomic annotation of a set of genes.

    If the annotation is not present locally, then mygene.info is used.
    All genome assemblies that are present in the latest version of
    Ensembl are supported by mygene.info.

    Parameters
    ----------
    genes : Iterable
        List of gene names or gene identifiers such as ensembl_id.
    genome : str
        Genome name

    Returns
    -------
    pandas.DataFrame with gene annotation.
    """

    # First try to find the genes in the local annotation installed by genomepy.
    gene_info = _local_gene_annotation(genes, genome)
    if gene_info is not None:
        return gene_info

    # Genes are not identified locally.
    # Retrieve the gene information using the mygene.info API
    logger.info(f"No local matching genes found for {genome}, trying mygene.info")

    # mygene.info only queries the most recent version of the Ensembl database
    # We can only safely continue if the local genome matched the Ensembl genome.
    # Even if the local genome was installed via Ensembl, we still need to check
    # if it is the same version
    result = ensembl_genome_info(genome)
    if result is None:
        return None

    # Run the actual query
    g = Genome(genome)
    logger.info("Querying mygene.info...")
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        genes,
        scopes="symbol,name,ensemblgene,entrezgene",
        fields="genomic_pos",
        species=g.tax_id,
        as_dataframe=True,
        verbose=False
    )

    if "notfound" in result and result.shape[1] == 1:
        logger.error("No matching genes found")
        sys.exit()

    if g.provider == "Ensembl":
        result = result.rename(columns={"genomic_pos.chr": "chrom"})
    else:
        # Ensembl, UCSC and NCBI chromosome names can all be different :-/
        logger.info("Local genome is not an Ensembl genome.")
        mapping = load_mapping(g.name)
        result = result.join(mapping, on="genomic_pos.chr")
        result = result.dropna(subset=["chrom"])

    # Convert genomic positions from string to integer
    result["genomic_pos.start"] = result["genomic_pos.start"].astype(int)
    result["genomic_pos.end"] = result["genomic_pos.end"].astype(int)

    # For each gene use the match with the highest score
    result = result.reset_index().sort_values("_score").groupby("query").last()

    # Map the Ensembl 1/-1 strand to +/- strand
    strand_df = pd.DataFrame(
        {"ens_strand": [-1, 1], "strand": ["-", "+"]}
    ).set_index("ens_strand")
    result = result.join(strand_df, on="genomic_pos.strand")

    # Select the correct columns and name them
    result = result.reset_index()[
        ["chrom", "genomic_pos.start", "genomic_pos.end", "query", "strand"]
    ]
    result.columns = [["chrom", "start", "end", "name", "strand"]]

    return result


if __name__ == "__main__":
    for genome in ["hg38", "Xenopus_tropicalis_v9.1"]:
        print(genome)
        print(gene_annotation(["TP53", "TP63", "GATA1", "GATA2"], genome))
