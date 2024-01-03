import os
import sys
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification
from lib.general_helpers.process_fastq import concatenate_fastq, extract_gz

# def the main function

def pipeline_directory(expedition_folder: str, consensus_method: str, windows: bool):
    input_base = 'assets/input'

    for root, dirs, files in os.walk(input_base):
        if root.endswith('fastq_pass'):
            i = 0
            for barcode_folder in dirs:
                print("Processing barcode folder: ", barcode_folder, " still in progress : ", i/len(dirs)*100, "%")
                barcode_path = os.path.join(root, barcode_folder)
                output_fastq = os.path.join(barcode_path, f"{barcode_folder}.fastq")
                intermediate_files = []

                # Extract .gz files and store paths for later deletion

                if not os.path.isfile(output_fastq):
                    for file in os.listdir(barcode_path):
                        if file.endswith('.gz'):
                            src_file = os.path.join(barcode_path, file)
                            dst_file = os.path.join(barcode_path, file[:-3])
                            extract_gz(src_file, dst_file)
                            intermediate_files.append(dst_file)

                    concatenate_fastq(barcode_path, output_fastq)

                    run_consensus(barcode_folder, output_fastq, consensus_method, barcode_path, wsl=windows)

                if os.path.getsize(os.path.join(barcode_path, f"{barcode_folder}_final_consensus.fasta")) != 0 :
                    print("Running identification for barcode folder: ", barcode_folder)
                    run_identification(barcode_folder, expedition_name=expedition_folder, input_path=barcode_path)


                # Delete intermediary files
                for file in intermediate_files:
                    os.remove(file)

                i += 1


    # identification


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pipeline_expedition.py <name_of_the_expedition> <name_of_the_consensus_method> <windows>")
        sys.exit(1)

    expedition_name = sys.argv[1]
    if len(sys.argv) == 4:
        name_of_the_consensus_method = sys.argv[2]
        windows = sys.argv[3]
    elif len(sys.argv) == 3:
        name_of_the_consensus_method = "80_20_best_sequence"
        windows = sys.argv[2]
    else:
        name_of_the_consensus_method = "80_20_best_sequence"
        windows = False
    

    pipeline_directory(expedition_name, name_of_the_consensus_method, windows)

if __name__ == "__main__":
    main()