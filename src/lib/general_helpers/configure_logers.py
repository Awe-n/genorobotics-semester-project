import os
import logging

def get_logger(name, output_dir, base_name):
    # Create a custom logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # Create handlers
    log_file = os.path.join(output_dir, f"{base_name}_{name}_pipeline_log.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    # Create formatters and add it to handlers
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    file_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)

    return logger

def configure_preprocessing_logger(output_dir, base_name):
    return get_logger("preprocessing", output_dir, base_name)

def configure_consensus_logger(output_dir, base_name):
    return get_logger("consensus", output_dir, base_name)

def configure_identification_logger(output_dir, base_name):
    return get_logger("identification", output_dir, base_name)