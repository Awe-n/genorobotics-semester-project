import os
from Bio import Entrez
from lib.general_helpers.run_command import run_command

def download_gene_sequences(gene_name: str, start: int = 0, end: int = 0, max_records: int = 1000000, batch_size: int = 10000):
    Entrez.email = "emilien.ordonneau@epfl.ch"
    query = f"{gene_name}[All Fields] AND (is_nuccore[filter] AND \"{start}\"[SLEN] : \"{end}\"[SLEN]))"
    print("Querying ", query)

    # First get count of total available records
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
    record_count = int(Entrez.read(handle)["Count"])
    handle.close()

    # Initialize variables
    fetched_records = 0
    all_sequences = ""

    # Fetch records in batches
    while fetched_records < min(record_count, max_records):
        print(f"Downloading from {fetched_records} to {fetched_records+batch_size} : {fetched_records/min(record_count,max_records)*100}%\r")
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

def make_blast_db(filename, db_name):
    command = f"makeblastdb -in {filename} -dbtype nucl -out {db_name}"
    run_command(command)

def check_blast_db(db_name):
    command = f"blastdbcmd -db {db_name} -info"
    return run_command(command)

def identification_pipeline_blastn(input_name: str, expedition_name: str = None, input_path: str = None, database: str = None, download: bool = False) -> list[tuple[str, str]]:
    databases = {"ITS", "matK", "psbA-trnH", "rbcL"} if database is None else {database}

    if download : 
        for db in databases:
            result, _ = check_blast_db(db)
            if not result.returncode == 0:
                if db == "matK":
                    filename = download_gene_sequences("matK", 750, 1500)
                elif db == "rbcL":
                    filename = download_gene_sequences("rbcL", 600, 1000)
                elif db == "psbA-trnH":
                    filename = download_gene_sequences("psbA-trnH", 400, 800)
                elif db == "ITS":
                    filename = download_gene_sequences("ITS")
                make_blast_db(filename, db)

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
    for db in databases:
        output_blastn_path = os.path.join(output_blastn, db + ".txt")
        blastn_cmd = f"blastn -query {input_path} -db {db} -out {output_blastn_path} -max_target_seqs 20 -outfmt 5"
        run_command(blastn_cmd)
        xml_files.append((output_blastn_path, db))

    print("You can find identification output at", output_blastn)

    return xml_files