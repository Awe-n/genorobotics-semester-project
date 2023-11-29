import logging
import time
import subprocess

def run_command(command):
    logging.info(f"Running command: {command}")
    start_time = time.time()
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    elapsed_time = time.time() - start_time
    logging.info(f"Command '{command}' took {elapsed_time:.2f} seconds to run.")
    if result.stdout:
        logging.info(result.stdout)
    if result.stderr:
        logging.error("Error: " + result.stderr)
    return result, elapsed_time