import os
import logging

def get_logger(name, output_dir, base_name):
    """
    Create and configure a custom logger.

    Args:
        name (str): The name of the logger.
        output_dir (str): The directory where the log file will be saved.
        base_name (str): The name of the input file the logger is associated with.

    Returns:
        logging.Logger: The configured logger instance.
    """
    # Create a custom logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # Prevent propagation to the root logger
    logger.propagate = False

    # Clear existing handlers (if any)
    if logger.hasHandlers():
        logger.handlers.clear()

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
    """
    Configure and return a logger for preprocessing.

    Args:
        output_dir (str): The directory where log files will be saved.
        base_name (str): The name of the input file the logger is associated with.

    Returns:
        logger: The configured logger for preprocessing.
    """
    return get_logger("preprocessing", output_dir, base_name)

def configure_consensus_logger(output_dir, base_name):
    """
    Configure and return a logger for consensus.

    Args:
        output_dir (str): The directory where the log files will be saved.
        base_name (str): The name of the input file the logger is associated with.
        
    Returns:
        logger: The configured consensus logger.
    """
    return get_logger("consensus", output_dir, base_name)

def configure_identification_logger(output_dir, base_name):
    """
    Configure and return a logger for identification.

    Args:
        output_dir (str): The directory where the log files will be saved.
        base_name (str): The name of the input file the logger is associated with.

    Returns:
        logger: The configured identification logger.
    """
    return get_logger("identification", output_dir, base_name)