from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_best_sequence import run_consensus_pipeline_80_20_best_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_80_20_longest_sequence import run_consensus_pipeline_80_20_longest_sequence
from lib.consensus.consensus_pipelines.consensus_pipeline_straightforward_best_sequence import run_consensus_pipeline_straightforward_best_sequence
from lib.general_helpers.configure_loggers import configure_consensus_logger
import os
import logging

def run_consensus(input_name: str, input_fastq_path: str, consensus_method: str, output_dir: str = None, logger : logging.Logger = None, windows: bool = False):
    """
    Runs the consensus pipeline based on the specified consensus method.

    Args:
        input_name (str): The name of the input.
        input_fastq_path (str): The path to the input FASTQ file.
        consensus_method (str): The method to use for generating the consensus.
        output_dir (str, optional): The output directory. If not provided, a default directory will be used.
        windows (bool, optional): Indicates whether the commands are to be ran in Windows Subsystem for Linux (WSL) environment.

    Raises:
        ValueError: If the specified consensus method is not recognized.
    """
    if output_dir == None : 
        output_dir = os.path.join("assets", "output", "consensus", input_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    if logger == None:
        logger = configure_consensus_logger(output_dir, input_name)
        print(f"Logging set up at {output_dir}/{input_name}_consensus_pipeline_log.log")

    logger.info("Running consensus pipeline...")
    
    # If statement to check which consensus method to run
    if consensus_method == "80_20_best_sequence":
        logger.info("Running consensus pipeline with 80_20_best_sequence method...")
        run_consensus_pipeline_80_20_best_sequence(input_name, input_fastq_path, output_dir, logger, windows)
    elif consensus_method == "80_20_longest_sequence":
        logger.info("Running consensus pipeline with 80_20_longest_sequence method...")
        run_consensus_pipeline_80_20_longest_sequence(input_name, input_fastq_path, output_dir, logger, windows)
    elif consensus_method == "straightforward_best_sequence":
        logger.info("Running consensus pipeline with streaming method...")
        run_consensus_pipeline_straightforward_best_sequence(input_name, input_fastq_path, output_dir, logger, windows)
    else:
        raise ValueError(f"Consensus method {consensus_method} not recognized.")