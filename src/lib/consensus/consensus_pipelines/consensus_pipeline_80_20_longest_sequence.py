import os
from lib.general_helpers.run_command import run_command
from lib.general_helpers.process_fastq import split_fastq
from lib.consensus.consensus_helpers.best_alignment import select_best_alignment
from Bio import SeqIO

def run_consensus_pipeline_80_20_longest_sequence(input_name: str, input_fastq_path: str, output_dir: str):
    input_fastq_filename = "rbcL_Qiagen_tomato_5000.fastq"
    input_fastq_path = f"assets/input/{input_fastq_filename}"
    base_name = os.path.splitext(input_fastq_filename)[0]

    output_base_dir = "assets/output"
    output_dir = os.path.join(output_base_dir, base_name)
    os.makedirs(output_dir, exist_ok=True)

    total_time_taken = 0

    # Split the file into top 20% and remaining 80%
    top_sequences_path, remaining_sequences_path = split_fastq(input_fastq_path, output_dir, base_name, percentile=20)

    # Step 1: Align and generate consensus from top 20% sequences
    top_paf_path = os.path.join(output_dir, f"{base_name}_top20_reads.paf")
    top_consensus_path = os.path.join(output_dir, f"{base_name}_top20_consensus.fasta")
    minimap2_command = f"minimap2 -x ava-ont {top_sequences_path} {top_sequences_path} > {top_paf_path}"
    print("Running read alignment with minimap2 on top 20% sequences...")
    _, minimap2_time = run_command(minimap2_command)
    total_time_taken += minimap2_time

    racon_command = f"racon -m 8 -x -6 -g -8 -w 500 {top_sequences_path} {top_paf_path} {top_sequences_path} > {top_consensus_path}"
    print("Generating consensus sequence with racon on top 20% sequences...")
    _, racon_time = run_command(racon_command)
    total_time_taken += racon_time

    # Step 1 bis: If the outputed racon file contains more than one sequence, we take the longest one
    consensus_sequences = list(SeqIO.parse(top_consensus_path, "fasta"))
    if len(consensus_sequences) > 1:
        print(f"Warning: {top_consensus_path} contains more than one sequence. Taking the longest one.")
        consensus_sequences.sort(key=lambda x: len(x), reverse=True)
        SeqIO.write(consensus_sequences[0], top_consensus_path, "fasta")

    # Step 2: Align the remaining 80% sequences to the top 20% consensus
    remaining_paf_path = os.path.join(output_dir, f"{base_name}_remaining80_reads.paf")
    final_consensus_path = os.path.join(output_dir, f"{base_name}_final_consensus.fasta")
    minimap2_command = f"minimap2 -x map-ont {top_consensus_path} {remaining_sequences_path} > {remaining_paf_path}"
    print("Running read alignment with minimap2 on remaining 80% sequences...")
    _, minimap2_time = run_command(minimap2_command)
    total_time_taken += minimap2_time

    # Step 3: Generate the final consensus sequence with racon
    racon_command = f"racon -m 8 -x -6 -g -8 -w 500 {remaining_sequences_path} {remaining_paf_path} {top_consensus_path} > {final_consensus_path}"
    print("Generating final consensus sequence with racon...")
    _, racon_time = run_command(racon_command)
    total_time_taken += racon_time

    # Print out the total time for each step
    print(f"Minimap2 alignment took {minimap2_time:.2f} seconds.")
    print(f"Total Racon iterations took {total_time_taken - minimap2_time:.2f} seconds.")

    # Print out the total time for the pipeline
    print(f"Total time taken for the pipeline: {total_time_taken:.2f} seconds.")

main()