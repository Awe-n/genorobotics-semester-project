from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_best_sequence import run_consensus_pipeline_80_20_best_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_longest_sequence import run_consensus_pipeline_80_20_longest_sequence

def run_consensus(input_name: str, input_fastq_path: str, output_dir: str, consensus_method: str = "80_20_best_sequence"):
    
    # If statement to check which consensus method to run
    if consensus_method == "80_20_best_sequence":
        run_consensus_pipeline_80_20_best_sequence(input_name, input_fastq_path, output_dir)
    elif consensus_method == "80_20_longest_sequence":
        run_consensus_pipeline_80_20_longest_sequence(input_name, input_fastq_path, output_dir)
    else:
        raise ValueError(f"Consensus method {consensus_method} not recognized.")