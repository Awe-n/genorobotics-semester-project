from lib.general_helpers.configure_loggers import configure_streaming_logger
from lib.streaming.streaming_pipelines.streaming_basic_pipeline import run_streaming_basic_pipeline
import os
import logging

def run_streaming(input_name: str, input_fastq_path: str, streaming_method: str, output_dir: str = None, logger : logging.Logger = None, wsl: bool = False):
    
    if output_dir == None : 
        output_dir = os.path.join("assets", "output", "streaming", input_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    if logger == None:
        logger = configure_streaming_logger(output_dir, input_name)
        print(f"Logging set up at {output_dir}/{input_name}_streaming_pipeline_log.log")

    logger.info("Running streaming pipeline...")

    # If statement to check which consensus method to run
    if streaming_method == "basic_streaming":
        logger.info("Running streaming pipeline with basic_streaming method...")
        run_streaming_basic_pipeline(input_name, input_fastq_path, output_dir, logger, wsl)
    else:
        raise ValueError(f"Streaming method {streaming_method} not recognized.")
