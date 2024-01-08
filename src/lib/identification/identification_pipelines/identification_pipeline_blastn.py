import os
from Bio import Entrez
from lib.general_helpers.run_command import run_command

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
    run_command(command, logger)

def check_blast_db(db_name, logger):
    """
    Check the information of a BLAST database.

    Args:
        db_name (str): The name of the BLAST database.

    Returns:
        str: The output of the 'blastdbcmd -db {db_name} -info' command.

    """
    command = f"blastdbcmd -db {db_name} -info"
    return run_command(command, logger)

def identification_pipeline_blastn(input_name: str, logger, expedition_name: str = None, input_path: str = None, database: str = None, download: bool = False) -> list[tuple[str, str]]:
    """
    Run the BLASTN identification pipeline.
    Works as follows : 
    1. Check if the database is available
    1b. If not, if download is True, download the gene sequences and make the database
    2. Run blastn for each specified database (ITS, matK, psbA-trnH, rbcL) (if not specified, run for all) and output the results to a XML file (assets/output/blastn/{input_name}/{database}.txt or assets/output/{expedition_name}/{input_name}/{database}.txt if expedition_name is not None)
    3. Return the path to the output XML file and the corresponding database name
    Note : as BLASTN exists in Windows and Linux, no need to specify wsl, in contrast to consensus

    Args:
        input_name (str): The name of the input.
        expedition_name (str, optional): The name of the expedition. Defaults to None.
        input_path (str, optional): The path to the input file. If None, use the default path (assets/output/consensus/{input_name}/{input_name}_final_consensus.fasta). Defaults to None.
        database (str, optional): The database to use for identification. If None, use all databases. Defaults to None.
        download (bool, optional): Whether to download the gene sequences if the database is not available. Defaults to False.

    Returns:
        list[tuple[str, str]]: A list of tuples containing the path to the output XML file and the corresponding database name.
    """
    databases = {"ITS", "matK", "psbA-trnH", "rbcL"} if database is None else {database}

    if download : 
        for curr_db in databases:
            result, _ = check_blast_db(curr_db, logger)
            if not result.returncode == 0:
                if curr_db == "matK":
                    filename = download_gene_sequences("matK", 750, 1500)
                elif curr_db == "rbcL":
                    filename = download_gene_sequences("rbcL", 600, 1000)
                elif curr_db == "psbA-trnH":
                    filename = download_gene_sequences("psbA-trnH", 400, 800)
                elif curr_db == "ITS":
                    filename = download_gene_sequences("ITS")
                make_blast_db(filename, curr_db, logger)

    if input_path is None :
        input_path = os.path.join("assets", "output", "consensus", input_name)

    input_path = os.path.join(input_path, f"{input_name}_final_consensus.fasta")

    output_blastn = os.path.join("assets", "output")
    if expedition_name is not None:
        output_blastn = os.path.join(output_blastn, expedition_name, input_name)
    else:
        output_blastn = os.path.join(output_blastn, "blastn", input_name)

    os.makedirs(output_blastn, exist_ok=True)

    xml_files = []
    for curr_db in databases:
        res = curr_db + ".txt"
        output_blastn_path = os.path.join(output_blastn, res)
        blastn_cmd = f"blastn -query {input_path} -db {curr_db} -out {output_blastn_path} -max_target_seqs 20 -outfmt 5"
        run_command(blastn_cmd, logger)
        xml_files.append((output_blastn_path, curr_db))

    return xml_files