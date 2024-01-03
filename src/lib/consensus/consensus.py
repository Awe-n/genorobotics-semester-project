from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_best_sequence import run_consensus_pipeline_80_20_best_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_longest_sequence import run_consensus_pipeline_80_20_longest_sequence
from lib.general_helpers.configure_logers import configure_consensus_logger
import os
import logging

def run_consensus(input_name: str, input_fastq_path: str, consensus_method: str, output_dir: str = None, wsl: bool = False):
    if output_dir == None : 
        output_dir = os.path.join("assets", "output", "consensus", input_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    log_file = configure_consensus_logger(output_dir, input_name)
    logging.info(f"Logging set up at {log_file}")

    logging.info("Running consensus pipeline...")
    print("Running consensus pipeline...")
    
    # If statement to check which consensus method to run
    if consensus_method == "80_20_best_sequence":
        logging.info("Running consensus pipeline with 80_20_best_sequence method...")
        print("Running consensus pipeline with 80_20_best_sequence method...")
        run_consensus_pipeline_80_20_best_sequence(input_name, input_fastq_path, output_dir, wsl)
    elif consensus_method == "80_20_longest_sequence":
        logging.info("Running consensus pipeline with 80_20_longest_sequence method...")
        print("Running consensus pipeline with 80_20_longest_sequence method...")
        run_consensus_pipeline_80_20_longest_sequence(input_name, input_fastq_path, output_dir, wsl)
    else:
        raise ValueError(f"Consensus method {consensus_method} not recognized.")
