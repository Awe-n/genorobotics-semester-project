import os
from Bio import Entrez
from lib.general_helpers.run_command import run_bash_command
from lib.identification.identification_helpers.make_blast_db import download_gene_sequences, make_blast_db

def check_blast_db(db_name, logger, windows: bool = False):
    """
    Check the information of a BLAST database.

    Args:
        db_name (str): The name of the BLAST database.

    Returns:
        str: The output of the 'blastdbcmd -db {db_name} -info' command.

    """
    command = f"blastdbcmd -db {db_name} -info"
    return run_bash_command(command, logger, windows)

def identification_pipeline_blastn(input_name: str, logger, expedition_name: str = None, input_path: str = None, output_dir: str = None,  database: str = None, download: bool = False, windows: bool = False) -> list[tuple[str, str]]:
    """
    Run the BLASTN identification pipeline.
    Works as follows : 
    1. Check if the database is available
    1b. If not, if download is True, download the gene sequences and make the database
    2. Run blastn for each specified database (ITS, matK, psbA-trnH, rbcL) (if not specified, run for all) and output the results to a XML file (assets/output/blastn/{input_name}/{database}.txt or assets/output/{expedition_name}/{input_name}/{database}.txt if expedition_name is not None)
    3. Return the path to the output XML file and the corresponding database name
    Note : as BLASTN exists in Windows and Linux, no need to specify windows, in contrast to consensus

    Args:
        input_name (str): The name of the input.
        expedition_name (str, optional): The name of the expedition. Defaults to None.
        input_path (str, optional): The path to the input file. If None, use the default path (assets/output/consensus/{input_name}/{input_name}_final_consensus.fasta). Defaults to None.
        database (str, optional): The database to use for identification. If None, use all databases. Defaults to None.
        download (bool, optional): Whether to download the gene sequences if the database is not available. Defaults to False.

    Returns:
        list[tuple[str, str]]: A list of tuples containing the path to the output XML file and the corresponding database name.
    """
    total_time_taken_blastn = 0

    databases = {"ITS", "matK", "psbA-trnH", "rbcL"} if database is None else {database}

    if download : 
        for curr_db in databases:
            result, _ = check_blast_db(curr_db, logger, windows)
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
        input_path = os.path.join("assets", "output", "post", input_name, "consensus")

    input_path = os.path.join(input_path, f"{input_name}_final_consensus.fasta")
    
    os.makedirs(output_dir, exist_ok=True)

    xml_files = []
    for curr_db in databases:
        res = curr_db + ".txt"
        output_blastn_path = os.path.join(output_dir, res)
        blastn_cmd = f"blastn -query {input_path} -db {curr_db} -out {output_blastn_path} -max_target_seqs 20 -outfmt 5"
        
        _, blastn_time = run_bash_command(blastn_cmd, logger, windows, True)
        logger.info(f"BLASTn command for database {curr_db} took {blastn_time:.2f} seconds.")

        total_time_taken_blastn += blastn_time

        xml_files.append((output_blastn_path, curr_db))

    return xml_files, total_time_taken_blastn