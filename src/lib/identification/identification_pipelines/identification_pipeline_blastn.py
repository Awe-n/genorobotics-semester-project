import os
from Bio import Entrez
from lib.general_helpers.run_command import run_command

def download_gene_sequences(gene_name: str, start: int= 0, end: int= 0):
    
    query = f"{gene_name}[Gene]" 
    if int != 0 : 
        query += f"AND {start}:{end}[Sequence Length]"

    search_handle = Entrez.esearch(db="nucleotide", term=query)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    id_list = search_results["IdList"]

    # Fetch sequences
    fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta")
    data = fetch_handle.read()
    fetch_handle.close()
    
    # Write to a file
    filename = f"{gene_name}_sequences.fasta"
    with open(filename, "w") as file:
        file.write(data)

    return filename

def make_blast_db(filename, db_name):
    command = f"makeblastdb -in {filename} -dbtype nucl -out {db_name}"
    run_command(command)

def check_blast_db(db_name):
    command = f"blastdbcmd -db {db_name} -info"
    return run_command(command)

def identification_pipeline_blastn(input_name: str, database: str = None):
    databases = {"ITS", "matK", "psbA-trnH", "rbcL"} if database is None else {database}

    for db in databases:
        if not check_blast_db(db):
            if db == "matK_seqs":
                filename = download_gene_sequences("matK", 750, 1500)
            elif db == "rbcL_seqs":
                filename = download_gene_sequences("rbcL", 600, 1000)
            elif db == "psbA-trnH_seqs":
                filename = download_gene_sequences("psbA-trnH", 400, 800)
            elif db == "ITS_seqs":
                filename = download_gene_sequences("ITS")
            make_blast_db(filename, db)

    input_fasta_path = os.path.join("assets", "output", input_name, input_name + "_final_consensus.fasta")
    output_blastn = os.path.join("assets", "output", "blastn", input_name)
    os.makedirs(output_blastn, exist_ok=True)


    for db in databases:
        output_blastn_path = os.path.join(output_blastn, db + ".txt")
        blastn_cmd = f"blastn -query {input_fasta_path} -db {db} -out {output_blastn_path} -max_target_seqs 5"
        run_command(blastn_cmd)

    print("You can find identification output at", output_blastn)