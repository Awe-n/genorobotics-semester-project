from lib.general_helpers.configure_logers import configure_identification_logger
from lib.identification.identification_pipelines.identification_pipeline_blastn import identification_pipeline_blastn
from lib.identification.identification_pipelines.identification_processing import get_best_species_from_xml
import os
import logging

def run_identification(input_name: str, expedition_name: str = None, input_path: str = None, output_dir: str = None, db: str = None):

    if (output_dir == None) : 
        output_dir = os.path.join("assets", "output", "blastn")

    # Configure logging
    logger = configure_identification_logger(output_dir, input_name)
    logger.info(f"Logging set up at {output_dir}/{input_name}_consensus_pipeline_log.log")

    logger.info("Running identification pipeline... /n")
    print("Running identification pipeline... /n")

    xml_files = identification_pipeline_blastn(logger, input_name, expedition_name, input_path, db)
    for xml_file, db in xml_files:
        best_species = get_best_species_from_xml(xml_file)
        logger.info(f"Best species for {db} is {best_species[0]} with alignment {best_species[1][0]} and evalue {best_species[1][1]}")
        print(f"Best species for {db} is {best_species[0]} with alignment {best_species[1][0]} and evalue {best_species[1][1]}")