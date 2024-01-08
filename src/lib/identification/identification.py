from lib.general_helpers.configure_loggers import configure_identification_logger
from lib.identification.identification_pipelines.identification_pipeline_blastn import identification_pipeline_blastn
from lib.identification.identification_pipelines.identification_processing import get_best_species_from_xml
import os
import logging

def write_results(output_dir: str, input_name: str, best_species_info: dict, dbs: set[str]):
    """
    Write the results of the identification pipeline to a file.

    Args:
        output_dir (str): The output directory.
        input_name (str): The name of the input.
        best_species_info (dict): The dictionary containing the best species information.
        dbs (set): The set of databases used.
    """

    with open(os.path.join(os.path.join(output_dir, input_name), f"{input_name}_identification_results.txt"), "w") as f:
        f.write(f"Input name: {input_name}\n")
        f.write(f"Databases used: {dbs}\n")
        f.write("\n")
        for db in dbs:
            f.write(f"Best species for {db} is {best_species_info[db]['species']} with alignment {best_species_info[db]['alignment']} and evalue {best_species_info[db]['evalue']}\n")

def run_identification(input_name: str, expedition_name: str = None, input_path: str = None, output_dir: str = None, database: str = None, logger: logging.Logger = None):
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
        output_dir = os.path.join("assets", "output", "blastn", input_name) if expedition_name == None else os.path.join("assets", "output", expedition_name)
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging
    if logger == None:
        logger = configure_identification_logger(output_dir, input_name)
        print(f"Logging set up at {output_dir}/{input_name}_identification_pipeline_log.log")

    logger.info("Running consensus pipeline... \n")

    best_species_info = {}

    xml_files = identification_pipeline_blastn(input_name, logger, expedition_name, input_path, database)
    logger.info(f"XML files : {xml_files}")
    for xml_file, db in xml_files:
        best_species = get_best_species_from_xml(xml_file)
        logger.info(f"Best species for {db} is {best_species[0]} with alignment {best_species[1][0]} and evalue {best_species[1][1]}")

        best_species_info[db] = {
            "species": best_species[0],
            "alignment": best_species[1][0],
            "evalue": best_species[1][1]
        }

    dbs = set(database) if database is not None else {"ITS", "matK", "psbA-trnH", "rbcL"}
    write_results(output_dir, input_name, best_species_info, dbs)

    return best_species_info