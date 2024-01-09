from lib.general_helpers.configure_loggers import configure_streaming_logger
from lib.streaming.streaming_pipelines.streaming_basic_pipeline import run_streaming_basic_pipeline
import os
import logging

def run_streaming(input_name: str, input_fastq_path: str, streaming_method: str, output_dir: str = None, logger : logging.Logger = None, db: str = None, wsl: bool = False, species_identification_percentage_dominance: float = 80.0, block_size: int = 500, minimum_block_amount_before_dominance_check: int = 5, consensus_method: str = "80_20_best_sequence", identification_method: str = "blastn"):
    """
    Main function to initiate and run the streaming pipeline based on the given parameters.
    :param input_name: Name of the input data.
    :param input_fastq_path: Path to the input FASTQ file.
    :param streaming_method: The method of streaming to be used.
    :param output_dir: Directory to save output data. Defaults to a structured path in 'assets/output/streaming'.
    :param logger: Logger object for logging the process. If None, a new logger is configured.
    :param db: Database name or path for identification.
    :param wsl: Boolean flag to indicate if Windows Subsystem for Linux (WSL) is used.
    :param species_identification_percentage_dominance: The threshold percentage for species dominance in identification.
    :param block_size: The number of reads in each block.
    :param minimum_block_amount_before_dominance_check: Minimum number of blocks to process before checking species dominance.
    :param consensus_method: Method used for consensus.
    :param identification_method: Method used for identification.
    """

    # Set up the output directory if not provided
    if output_dir is None: 
        output_dir = os.path.join("assets", "output", "streaming", input_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    if logger is None:
        logger = configure_streaming_logger(output_dir, input_name)
        print(f"Logging set up at {output_dir}/{input_name}_streaming_pipeline_log.log")

    # Log the details of the streaming pipeline
    logger.info("Details of the streaming pipeline:" + "\n")
    logger.info("Streaming method: " + streaming_method + "\n")
    logger.info("The block size is the number of reads in each block.")
    logger.info("Streaming block size: " + str(block_size) + "\n")
    logger.info("Species identification percentage dominance: " + str(species_identification_percentage_dominance * 100) + "%")
    logger.info("The species identification percentage dominance is the minimum percentage of blocks that must identify the same species for the species to be considered as the consensus species.")
    logger.info("When the species identification percentage dominance is reached, the streaming pipeline will stop." + "\n")
    logger.info(f"--------------------------------------")
    logger.info("Running streaming pipeline...")

    # Execute the appropriate streaming pipeline based on the specified method
    if streaming_method == "basic_streaming":
        logger.info("Running streaming pipeline with basic_streaming method...")
        run_streaming_basic_pipeline(input_name, input_fastq_path, output_dir, logger, db, wsl, species_identification_percentage_dominance, block_size, minimum_block_amount_before_dominance_check, consensus_method, identification_method)
    else:
        raise ValueError(f"Streaming method {streaming_method} not recognized.")
