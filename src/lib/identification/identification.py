from lib.general_helpers.configure_loggers import configure_identification_logger
from lib.identification.identification_pipelines.identification_pipeline_blastn import identification_pipeline_blastn
from lib.identification.identification_pipelines.identification_processing import get_best_species_from_xml
import os
import logging

def run_identification(input_name: str, expedition_name: str = None, input_path: str = None, output_dir: str = None, db: str = None, logger: logging.Logger = None, identification_method: str = "blastn"):
    """
    Run the identification pipeline by following these steps:
    1. Run the BLASTN identification pipeline.
    2. Get the best species from the resulting XML file.

    Args:
        input_name (str): The name of the input.
        expedition_name (str, optional): The name of the expedition. Defaults to None.
        input_path (str, optional): The path of the input. Defaults to None.
        output_dir (str, optional): The output directory. Defaults to None.
        db (str, optional): The database to be used by BLASTN. Defaults to None.
    """

    if (output_dir == None) : 
        output_dir = os.path.join("assets", "output", "post", input_name, "identification") if expedition_name == None else os.path.join("assets", "output", expedition_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    if logger == None:
        logger = configure_identification_logger(output_dir, input_name)
        print(f"Logging set up at {output_dir}/{input_name}_identification_pipeline_log.log")

    logger.info("Running consensus pipeline... \n")

    # If statement to check which identification method to run
    if identification_method == "blastn":
        logger.info("Running identification pipeline with blastn method... \n")

        best_species_info = {}

        xml_files, total_time_taken_blastn = identification_pipeline_blastn(input_name, logger, expedition_name, input_path, output_dir, db)

        logger.info(f"XML files : {xml_files}")
        for xml_file, db in xml_files:
            best_species = get_best_species_from_xml(xml_file)
            logger.info(f"Best species for {db} is {best_species[0]} with alignment {best_species[1][0]} and evalue {best_species[1][1]}")

            best_species_info[db] = {
                "species": best_species[0],
                "alignment": best_species[1][0],
                "evalue": best_species[1][1]
            }

        return best_species_info, total_time_taken_blastn
    else:
        raise ValueError(f"Identification method {identification_method} not recognized.")