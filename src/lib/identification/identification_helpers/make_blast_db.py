from Bio import Entrez
from lib.general_helpers.run_command import run_bash_command

def download_gene_sequences(gene_name: str, logger, start: int = 0, end: int = 0, max_records: int = 1000000, batch_size: int = 10000):
    """
    Downloads gene sequences from the NCBI nucleotide database based on the provided gene name, start and end positions,
    maximum number of records, and batch size.

    Args:
        gene_name (str): The gene name to search for.
        start (int, optional): The start position of the gene sequence. Defaults to 0.
        end (int, optional): The end position of the gene sequence. Defaults to 0.
        max_records (int, optional): The maximum number of records to download. Defaults to 1000000.
        batch_size (int, optional): The batch size for downloading records. Defaults to 10000.

    Returns:
        str: The filename of the downloaded gene sequences.
    """
    # Set email for Entrez
    Entrez.email = None
    query = f"{gene_name}[All Fields] AND (is_nuccore[filter] AND \"{start}\"[SLEN] : \"{end}\"[SLEN]))"
    logger.info("Querying ", query)

    # First get count of total available records
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
    record_count = int(Entrez.read(handle)["Count"])
    handle.close()

    # Initialize variables
    fetched_records = 0
    all_sequences = ""

    # Fetch records in batches
    while fetched_records < min(record_count, max_records):
        logger.info(f"Downloading from {fetched_records} to {fetched_records+batch_size} : {fetched_records/min(record_count,max_records)*100}%\r")
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=batch_size, retstart=fetched_records)
        search_results = Entrez.read(handle)
        handle.close()
        id_list = search_results["IdList"]

        handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        data = handle.read()
        handle.close()

        all_sequences += data
        fetched_records += len(id_list)

    # Write to a file
    filename = f"{gene_name}_sequences.fasta"
    with open(filename, "w") as file:
        file.write(all_sequences)

    return filename

def make_blast_db(filename, db_name, logger):
    """
    Create a BLAST database using the specified input file.

    Args:
        filename (str): The path to the input file.
        db_name (str): The name of the BLAST database to be created.

    Returns:
        None
    """
    command = f"makeblastdb -in {filename} -dbtype nucl -out {db_name}"
    run_bash_command(command, logger)