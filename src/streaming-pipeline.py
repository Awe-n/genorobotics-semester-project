import os
import sys
from lib.streaming.streaming import run_streaming

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_fastq_filename> <windows (True/False)> <db> <streaming_method> <consensus_method> <identification_method>")
        sys.exit(1)

    # Assign command-line arguments to variables
    input_fastq_filename = sys.argv[1]
    windows = sys.argv[2].lower() in ["true", "t", "yes", "y"]
    db = sys.argv[3]
    streaming_method = sys.argv[4]
    consensus_method = sys.argv[5]
    identification_method = sys.argv[6]

    # Default values for other parameters (can be modified as needed)
    species_identification_percentage_dominance = 0.7
    block_size = 250
    minimum_block_amount_before_dominance_check = 5

    input_fastq_path = os.path.join("assets", "input", input_fastq_filename)
    base_name = os.path.splitext(input_fastq_filename)[0]

    run_streaming(input_name=base_name,
                  input_fastq_path=input_fastq_path,
                  streaming_method=streaming_method,
                  db=db,
                  windows=windows,
                  species_identification_percentage_dominance=species_identification_percentage_dominance,
                  block_size=block_size,
                  minimum_block_amount_before_dominance_check=minimum_block_amount_before_dominance_check,
                  consensus_method=consensus_method,
                  identification_method=identification_method)

if __name__ == "__main__":
    main()