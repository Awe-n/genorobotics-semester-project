import os
import sys
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification
from lib.general_helpers.process_fastq import concatenate_fastq, extract_gz
from lib.general_helpers.configure_loggers import get_logger

# def the main function

def pipeline_directory(expedition_folder: str, consensus_method: str, windows: bool):
    input_base = os.path.join('assets/input', expedition_folder)
    log_path = os.path.join('assets/log/', expedition_folder)
    os.makedirs(log_path, exist_ok=True)
    logger = get_logger(expedition_folder, log_path, expedition_folder)
    print(f"Logging set up at {log_path}/{expedition_folder}_pipeline_log.log")

    for root, dirs, files in os.walk(input_base):
        if root.endswith('fastq_pass'):
            i = 0
            for barcode_folder in dirs:
                # Print the current progress of the pipeline
                logger.info(f"Processing barcode folder: {barcode_folder}, still in progress : {i/len(dirs)*100:.2f}%")
                print(end='\x1b[2K')
                print(f"\rProcessing barcode folder: {barcode_folder}, still in progress : {i/len(dirs)*100:.2f}%", end='', flush=True)
                barcode_path = os.path.join(root, barcode_folder)
                output_fastq = os.path.join(barcode_path, f"{barcode_folder}.fastq")
                intermediate_files = []

                # Extract .gz files and concatenate them into a single .fastq file, then generate the consensus sequence from it if it doesn't exist yet
                if not os.path.isfile(output_fastq):
                    logger.info(f"Extracting and concatenating files for barcode folder {barcode_folder}")
                    for file in os.listdir(barcode_path):
                        if file.endswith('.gz'):
                            src_file = os.path.join(barcode_path, file)
                            dst_file = os.path.join(barcode_path, file[:-3])
                            extract_gz(src_file, dst_file)
                            intermediate_files.append(dst_file)

                    concatenate_fastq(barcode_path, output_fastq)

                    logger.info(f"Running consensus for barcode folder {barcode_folder}")
                    run_consensus(barcode_folder, output_fastq, consensus_method, barcode_path, wsl=windows)
                else :
                    logger.info(f"Consensus already computed for barcode folder {barcode_folder}")

                # Run identification if the consensus sequence exists
                if os.path.getsize(os.path.join(barcode_path, f"{barcode_folder}_final_consensus.fasta")) != 0 :
                    logger.info(f"Running identification for barcode folder {barcode_folder}")
                    run_identification(barcode_folder, expedition_name=expedition_folder, input_path=barcode_path, logger=logger)
                else :
                    logger.info(f"No consensus found for barcode folder {barcode_folder}")


                # Delete intermediary files
                for file in intermediate_files:
                    os.remove(file)

                i += 1
    print(end='\x1b[2K')
    print("\rProcessing barcode folder: Done! 100.00%")


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pipeline_expedition.py <name_of_the_expedition> <windows> <name_of_the_consensus_method>")
        sys.exit(1)

    expedition_name = sys.argv[1]
    windows = False
    name_of_the_consensus_method = "80_20_best_sequence"

    if len(sys.argv) >= 3:
        windows = True if sys.argv[2].lower() in ["true", "t", "yes", "y"] else False
    if len(sys.argv) >= 4:
        name_of_the_consensus_method = sys.argv[2]
    

    pipeline_directory(expedition_name, name_of_the_consensus_method, windows)

if __name__ == "__main__":
    main()