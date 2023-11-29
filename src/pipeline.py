import os
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification

# def the main function

def main():
    # initialize variables
    input_fastq_filename = "rbcL_Qiagen_tomato_5000.fastq"
    input_fastq_path = f"assets/input/{input_fastq_filename}"
    base_name = os.path.splitext(input_fastq_filename)[0]

    output_base_dir = "assets/output"
    output_dir = os.path.join(output_base_dir, base_name)
    os.makedirs(output_dir, exist_ok=True)

    # preprocessing

    # run_preprocessing(input_fastq_filename, input_fastq_path, output_dir)

    # consensus

    consensus_method = "80_20_best_sequence"
    run_consensus(base_name, input_fastq_path, output_dir, consensus_method)

    # identification

    run_identification(base_name, output_dir)


if __name__ == "__main__":
    main()