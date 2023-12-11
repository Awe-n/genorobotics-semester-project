from lib.identification.identification_pipelines.identification_pipeline_blastn import identification_pipeline_blastn
from lib.general_helpers.configure_logers import configure_identification_logger
import os
import logging

def run_identification(input_name: str, output_dir: str = None, db: str = None):
    if (output_dir == None) : 
        output_dir = os.path.join("assets", "output", "blastn")

    # Configure logging
    log_file = configure_identification_logger(output_dir, input_name)
    logging.info(f"Logging set up at {log_file}")

    logging.info("Running consensus pipeline... \n")
    print("Running consensus pipeline... \n")

    identification_pipeline_blastn(input_name, db)