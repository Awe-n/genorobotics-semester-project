import os
from datetime import datetime

def write_identification_results(output_dir: str, input_name: str, best_species_info: dict, dbs: set[str]):
    """
    Write the results of the identification pipeline to a file.

    Args:
        output_dir (str): The output directory.
        input_name (str): The name of the input.
        best_species_info (dict): The dictionary containing the best species information.
        dbs (set): The set of databases used.
    """

    if not isinstance(dbs, set):
        raise ValueError("dbs should be a set of database names.")

    # Format the current date and time for the log filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    with open(os.path.join(output_dir, f"{input_name}_identification_results_{timestamp}.txt"), "w") as f:
        f.write(f"Input name: {input_name}\n")
        f.write(f"Databases used: {', '.join(dbs)}\n")
        f.write("\n")

        for db in dbs:
            if db in best_species_info:
                f.write(f"Best species for {db} is {best_species_info[db]['species']} with alignment {best_species_info[db]['alignment']} and evalue {best_species_info[db]['evalue']}\n")
            else:
                f.write(f"No information available for database {db}\n")
