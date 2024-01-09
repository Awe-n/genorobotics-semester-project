import os
from lib.general_helpers.run_command import run_command
from lib.consensus.consensus_helpers.best_alignment import select_best_alignment
from Bio import SeqIO

def run_consensus_pipeline_straightforward_best_sequence(input_name: str, input_fastq_path: str, output_dir: str, logger, windows: bool = False):
    """
    Runs a straightforward consensus pipeline with the following steps:
    1. Aligns all reads using minimap2.
    2. Generates a consensus sequence with racon.
    3. If the output racon file contains more than one sequence, selects the best alignment.

    Args:
        input_name (str): The name of the input file.
        input_fastq_path (str): The path to the input FASTQ file.
        output_dir (str): The directory to save the output files.
        logger: The logger object for logging messages.
        windows (bool, optional): Indicates whether running on Windows Subsystem for Linux (windows). Defaults to False.
    """
    total_time_taken = 0
    total_time_taken_minimap2 = 0
    total_time_taken_racon = 0

    # Step 1: Align all reads using minimap2
    paf_path = os.path.join(output_dir, f"{input_name}_reads.paf")
    minimap2_command = f"minimap2 -x ava-ont {input_fastq_path} {input_fastq_path} > {paf_path}"
    logger.info("Running read alignment with minimap2...")
    _, minimap2_time = run_command(minimap2_command, logger, windows)
    total_time_taken_minimap2 += minimap2_time

    # Step 2: Generate consensus sequence with racon
    consensus_path = os.path.join(output_dir, f"{input_name}_final_consensus.fasta")
    racon_command = f"racon -m 8 -x -6 -g -8 -w 500 {input_fastq_path} {paf_path} {input_fastq_path} > {consensus_path}"
    logger.info("Generating consensus sequence with racon...")
    _, racon_time = run_command(racon_command, logger, windows)
    total_time_taken_racon += racon_time

    # Step 3: Select the best alignment if multiple sequences are present
    consensus_sequences = list(SeqIO.parse(consensus_path, "fasta"))
    if len(consensus_sequences) > 1:
        logger.info(f"Multiple sequences found in {consensus_path}. Selecting the best alignment...")
        best_sequence = select_best_alignment(consensus_sequences, consensus_sequences[0])  # Choose the first sequence as the reference
        SeqIO.write(best_sequence, consensus_path, "fasta")

    # Delete intermediate files
    logger.info("Deleting intermediate files...")
    os.remove(paf_path)

    # Log the total time for each step
    logger.info(f"Minimap2 alignment took {total_time_taken_minimap2:.2f} seconds.")
    logger.info(f"Racon consensus generation took {total_time_taken_racon:.2f} seconds.")

    # Log the total time for the pipeline
    total_time_taken = total_time_taken_minimap2 + total_time_taken_racon
    logger.info(f"Total time taken for the consensus pipeline: {total_time_taken:.2f} seconds.")

    return total_time_taken, total_time_taken_minimap2, total_time_taken_racon