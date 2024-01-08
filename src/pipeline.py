import os
import sys
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification

# def the main function

def pipeline_single(input_fastq_filename: str, db: str, windows: bool):

    input_fastq_path = f"assets/input/{input_fastq_filename}"
    base_name = os.path.splitext(input_fastq_filename)[0]

    # preprocessing

    # run_preprocessing(input_fastq_filename, input_fastq_path, output_dir)

    # consensus

    consensus_method = "80_20_best_sequence"
    run_consensus(base_name, input_fastq_path, consensus_method, windows=windows)

    # identification

    run_identification(base_name, db=db)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 pipeline.py <name_of_the_input_fastq_file> <windows (True/False, Optional)> <db (Optional)>")
        sys.exit(1)

    input_fastq_filename = sys.argv[1]
    windows = False
    db = None
    if len(sys.argv) >= 3:
        windows = True if sys.argv[2].lower() in ["true", "t", "yes", "y"] else False
    if len(sys.argv) >= 4:
        db = sys.argv[3]

    pipeline_single(input_fastq_filename, db, windows)


if __name__ == "__main__":
    main()