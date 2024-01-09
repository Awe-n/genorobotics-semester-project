import subprocess
import time

def run_bash_command(command: str, logger, windows: bool = False):
    """
    Runs a command in the shell and logs the command, execution time, stdout, and stderr.

    Args:
        command (str): The command to be executed.
        logger: The logger object used for logging.
        windows (bool, optional): Specifies whether the environment is Windows (True) or UNIX (False). 
                                  Defaults to False.

    Returns:
        tuple: A tuple containing the result object and the elapsed time in seconds.
    """
    if windows:
        # Adjust the command for Windows environment
        full_command = ("windows " + command).replace('\\', '/')
    else:
        # Use "bash -l -c" to run the command within the context of the bash profile on UNIX systems
        full_command = f'bash -l -c "{command}"'

    logger.info(f"Running command: {full_command}")
    start_time = time.time()
    result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    elapsed_time = time.time() - start_time
    logger.info(f"Command '{full_command}' took {elapsed_time:.2f} seconds to run.")
    if result.stdout:
        logger.info(result.stdout)
    if result.stderr:
        logger.error("Error: " + result.stderr)
    return result, elapsed_time

def run_command(command: str, logger, windows: bool = False):
    """
    Runs a command in the shell and logs the command, execution time, stdout, and stderr.

    Args:
        command (str): The command to be executed.
        logger: The logger object used for logging.
        windows (bool, optional): Specifies whether to run the command in Windows Subsystem for Linux (windows). 
                              Defaults to False.

    Returns:
        tuple: A tuple containing the result object and the elapsed time in seconds.
    """
    if windows : 
        command = ("windows " + command).replace('\\','/')
    logger.info(f"Running command: {command}")
    start_time = time.time()
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    elapsed_time = time.time() - start_time
    logger.info(f"Command '{command}' took {elapsed_time:.2f} seconds to run.")
    if result.stdout:
        logger.info(result.stdout)
    if result.stderr:
        logger.error("Error: " + result.stderr)
    return result, elapsed_time