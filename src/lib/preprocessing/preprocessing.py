from lib.general_helpers.configure_logers import configure_preprocessing_logger
import logging

def run_prepocessing(input_fastq_filename, input_fastq_path, output_dir):

    # Configure logging
    log_file = configure_preprocessing_logger(output_dir, input_fastq_filename)
    logging.info(f"Logging set up at {log_file}")

    logging.info("Running preprocessing pipeline... /n")
    print("Running preprocessing pipeline... /n")

    # Run preprocessing pipeline
    # preprocessing_pipeline(input_fastq_filename, input_fastq_path, output_dir)