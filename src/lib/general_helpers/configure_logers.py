import os
import logging

def configure_preprocessing_logger(output_dir, base_name):
    log_file = os.path.join(output_dir, f"{base_name}_preprocessing_pipeline_log.log")
    logging.basicConfig(filename=log_file, level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    return log_file

def configure_consensus_logger(output_dir, base_name):
    log_file = os.path.join(output_dir, f"{base_name}_consensus_pipeline_log.log")
    logging.basicConfig(filename=log_file, level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    return log_file

def configure_identification_logger(output_dir, base_name):
    log_file = os.path.join(output_dir, f"{base_name}_identification_pipeline_log.log")
    logging.basicConfig(filename=log_file, level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    return log_file