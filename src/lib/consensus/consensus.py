from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_best_sequence import run_consensus_pipeline_80_20_best_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_longest_sequence import run_consensus_pipeline_80_20_longest_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_streaming import run_consensus_streaming_pipeline
from lib.consensus.consensus_pipelines.consensus_pipeline_streaming_wip import run_consensus_streaming_pipeline_wip
from lib.general_helpers.configure_logers import configure_consensus_logger
import os
import logging

def run_consensus(input_name: str, input_fastq_path: str, consensus_method: str, output_dir: str = None, wsl: bool = False):
    if output_dir == None : 
        output_dir = os.path.join("assets", "output", "consensus", input_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    logger = configure_consensus_logger(output_dir, input_name)
    logger.info(f"Logging set up at {output_dir}/{input_name}_consensus_pipeline_log.log")

    logger.info("Running consensus pipeline... /n")
    print("Running consensus pipeline... /n")
    
    # If statement to check which consensus method to run
    if consensus_method == "80_20_best_sequence":
        logger.info("Running consensus pipeline with 80_20_best_sequence method...")
        print("Running consensus pipeline with 80_20_best_sequence method...")
        run_consensus_pipeline_80_20_best_sequence(input_name, input_fastq_path, output_dir, logger, wsl)
    elif consensus_method == "80_20_longest_sequence":
        logger.info("Running consensus pipeline with 80_20_longest_sequence method...")
        print("Running consensus pipeline with 80_20_longest_sequence method...")
        run_consensus_pipeline_80_20_longest_sequence(input_name, input_fastq_path, output_dir, logger, wsl)
    elif consensus_method == "streaming":
        logger.info("Running consensus pipeline with streaming method...")
        print("Running consensus pipeline with streaming method...")
        run_consensus_streaming_pipeline(input_name, input_fastq_path, output_dir, logger, wsl)
    elif consensus_method == "streaming_wip":
        logger.info("Running consensus pipeline with streaming method WORK IN PROGRESS...")
        print("Running consensus pipeline with streaming method WORK IN PROGRESS...")
        run_consensus_streaming_pipeline_wip(input_name, input_fastq_path, output_dir, logger, wsl)
    else:
        raise ValueError(f"Consensus method {consensus_method} not recognized.")
