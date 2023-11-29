import os
from ...general_helpers.run_command import run_command

def identification_pipeline_blastn(input_name: str, database: str = None):
    input_fasta_path = os.path.join("assets", "output", input_name, input_name + "_final_consensus.fasta")
    output_blastn = os.path.join("assets", "output", "blastn", input_name)
    os.makedirs(output_blastn, exist_ok=True)

    databases = {"ITS_seqs", "matK_seqs", "psbA-trnH_seqs", "rbcL_seqs"} if database is None else {database}

    for db in databases:
        output_blastn_path = os.path.join(output_blastn, db + ".txt")
        blastn_cmd = f"blastn -query {input_fasta_path} -db {db} -out {output_blastn_path} -max_target_seqs 5"
        run_command(blastn_cmd)

    print("You can find identification output at", output_blastn)