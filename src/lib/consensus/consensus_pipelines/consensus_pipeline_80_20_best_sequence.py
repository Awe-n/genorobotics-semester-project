import os
from lib.general_helpers.run_command import run_command
from lib.general_helpers.process_fastq import split_fastq
from lib.consensus.consensus_helpers.best_alignment import select_best_alignment
from Bio import SeqIO

def run_consensus_pipeline_80_20_best_sequence(input_name: str, input_fastq_path: str, output_dir: str, logger, windows: bool = False):
    """
    Runs the consensus pipeline with the following steps:
    1. Splits the input FASTQ file into top 20% and remaining 80% sequences.
    2. Aligns and generates consensus from the top 20% sequences using minimap2 and racon.
    3. If the output racon file contains more than one sequence, selects the best alignment.
    4. Aligns the remaining 80% sequences to the top 20% consensus using minimap2.
    5. Generates the final consensus sequence with racon.

    Args:
        input_name (str): The name of the input file.
        input_fastq_path (str): The path to the input FASTQ file.
        output_dir (str): The directory to save the output files.
        logger: The logger object for logging messages.
        windows (bool, optional): Indicates whether running on Windows. Defaults to False.
    """
    total_time_taken = 0
    total_time_taken_minimap2 = 0
    total_time_taken_racon = 0

    # Split the file into top 20% and remaining 80%
    top_sequences_path, remaining_sequences_path = split_fastq(input_fastq_path, output_dir, input_name, percentile=20)

    # Step 1: Align and generate consensus from top 20% sequences
    top_paf_path = os.path.join(output_dir, f"{input_name}_top20_reads.paf")
    top_consensus_path = os.path.join(output_dir, f"{input_name}_top20_consensus.fasta")
    minimap2_command = f"minimap2 -x ava-ont {top_sequences_path} {top_sequences_path} > {top_paf_path}"
    logger.info("Running read alignment with minimap2 on top 20% sequences...")
    _, minimap2_time = run_command(minimap2_command, logger, windows)
    total_time_taken_minimap2 += minimap2_time

    racon_command = f"racon -m 8 -x -6 -g -8 -w 500 {top_sequences_path} {top_paf_path} {top_sequences_path} > {top_consensus_path}"
    logger.info("Generating consensus sequence with racon on top 20% sequences...")
    _, racon_time = run_command(racon_command, logger, windows)
    total_time_taken_racon += racon_time

    # Step 1 bis: If the output racon file contains more than one sequence, select the best alignment
    consensus_sequences = list(SeqIO.parse(top_consensus_path, "fasta"))
    if len(consensus_sequences) > 1:
        logger.info(f"Multiple sequences found in {top_consensus_path}. Selecting the best alignment...")
        best_sequence = select_best_alignment(consensus_sequences, consensus_sequences[0])  # Choose the first sequence as the reference
        SeqIO.write(best_sequence, top_consensus_path, "fasta")

    # Step 2: Align the remaining 80% sequences to the top 20% consensus
    remaining_paf_path = os.path.join(output_dir, f"{input_name}_remaining80_reads.paf")
    final_consensus_path = os.path.join(output_dir, f"{input_name}_final_consensus.fasta")
    minimap2_command = f"minimap2 -x map-ont {top_consensus_path} {remaining_sequences_path} > {remaining_paf_path}"
    logger.info("Running read alignment with minimap2 on remaining 80% sequences...")
    _, minimap2_time = run_command(minimap2_command, logger, windows)
    total_time_taken_minimap2 += minimap2_time

    # Step 3: Generate the final consensus sequence with racon
    racon_command = f"racon -m 8 -x -6 -g -8 -w 500 {remaining_sequences_path} {remaining_paf_path} {top_consensus_path} > {final_consensus_path}"
    logger.info("Generating final consensus sequence with racon...")
    _, racon_time = run_command(racon_command, logger, windows)
    total_time_taken_racon += racon_time    

    # Delete intermediate files
    logger.info("Deleting intermediate files...")
    os.remove(top_paf_path)
    os.remove(remaining_paf_path)
    os.remove(top_sequences_path)
    os.remove(remaining_sequences_path)

    # Log the total time for each step
    logger.info(f"Minimap2 alignment took {total_time_taken_minimap2:.2f} seconds.")
    logger.info(f"Total Racon iterations took {total_time_taken_racon:.2f} seconds.")

    # Log the total time for the pipeline
    total_time_taken = total_time_taken_minimap2 + total_time_taken_racon
    logger.info(f"Total time taken for the consensus pipeline: {total_time_taken:.2f} seconds.")
