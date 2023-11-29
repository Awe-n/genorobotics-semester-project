from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_best_sequence import run_consensus_pipeline_80_20_best_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_longest_sequence import run_consensus_pipeline_80_20_longest_sequence
import os

def run_consensus(input_name: str, input_fastq_path: str, output_dir: str, consensus_method: str):

    os.makedirs(output_dir, exist_ok=True)

    print("Running consensus pipeline...")
    
    # If statement to check which consensus method to run
    if consensus_method == "80_20_best_sequence":
        print("Running consensus pipeline with 80_20_best_sequence method...")
        run_consensus_pipeline_80_20_best_sequence(input_name, input_fastq_path, output_dir)
    elif consensus_method == "80_20_longest_sequence":
        print("Running consensus pipeline with 80_20_longest_sequence method...")
        run_consensus_pipeline_80_20_longest_sequence(input_name, input_fastq_path, output_dir)
    else:
        raise ValueError(f"Consensus method {consensus_method} not recognized.")