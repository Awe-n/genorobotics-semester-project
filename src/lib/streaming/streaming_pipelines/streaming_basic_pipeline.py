from lib.streaming.streaming_helpers.random_sampling import random_sampling_blocks, save_block_as_fastq, read_streaming_fastq
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification
from lib.general_helpers.configure_loggers import configure_consensus_logger
from lib.general_helpers.configure_loggers import configure_identification_logger
from lib.identification.identification_helpers.write_identification_results import write_identification_results
import os

def setup_and_save_block(block, block_index, output_dir, streaming_logger):
    # Creates a directory for each block and saves the block as a FASTQ file.
    block_output_dir = os.path.join(output_dir, f"block_{block_index}")
    os.makedirs(block_output_dir, exist_ok=True)
    save_block_as_fastq(block, block_index, block_output_dir)
    streaming_logger.info(f"Block {block_index} saved with {len(block)} reads")
    return block_output_dir

def run_consensus_for_block(block_input_name, block_fastq_path, consensus_method, block_output_dir, wsl, streaming_logger):
    # Creates a consensus directory within the block directory and runs the consensus pipeline.
    consensus_block_output_dir = os.path.join(block_output_dir, "consensus")
    os.makedirs(consensus_block_output_dir, exist_ok=True)
    consensus_logger = configure_consensus_logger(consensus_block_output_dir, block_input_name)
    streaming_logger.info(f"Running consensus pipeline for block {block_input_name}...")
    return run_consensus(block_input_name, block_fastq_path, consensus_method, consensus_block_output_dir, consensus_logger, wsl)

def run_identification_for_block(block_input_name, block_output_dir, db, identification_method, streaming_logger, consensus_block_output_dir, windows: bool = False):
    # Creates an identification directory within the block directory and runs the identification pipeline.
    identification_block_output_dir = os.path.join(block_output_dir, "identification")
    os.makedirs(identification_block_output_dir, exist_ok=True)
    identification_logger = configure_identification_logger(identification_block_output_dir, block_input_name)
    streaming_logger.info(f"Running identification pipeline for block {block_input_name}...")
    return run_identification(input_name=block_input_name, input_path=consensus_block_output_dir, output_dir=identification_block_output_dir, db=db, logger=identification_logger, identification_method=identification_method, windows=windows)

def check_species_dominance(blastn_result_list, species_identification_percentage_dominance, streaming_logger, output_dir, input_name, db=None):
    """
    Checks for the dominance of a species in the identification results. Logs details about the dominant species
    and other identified species, including their average e-values.
    """
    species_list = [result['species'] for blastn_result in blastn_result_list for db, result in blastn_result.items()]
    evalues = {species: [] for species in set(species_list)}  # Store e-values for each species
    species_count = {species: 0 for species in set(species_list)}

    # Count species occurrences and store e-values
    for blastn_result in blastn_result_list:
        for db, result in blastn_result.items():
            species = result['species']
            species_count[species] += 1
            evalues[species].append(result['evalue'])

    total_identifications = len(species_list)
    dominant_species = None
    dominant_percentage = 0

    # Find dominant species and calculate average e-value for each species
    for species, count in species_count.items():
        percentage = (count / total_identifications) * 100
        avg_evalue = sum(evalues[species]) / len(evalues[species])
        avg_alignment = sum([result['alignment'] for blastn_result in blastn_result_list for db, result in blastn_result.items() if result['species'] == species]) / count
        if percentage > dominant_percentage:
            dominant_percentage = percentage
            dominant_species = species
        streaming_logger.info(f"Species {species} identified {percentage:.2f}% of times with average alignment {avg_alignment:.2f} and average e-value {avg_evalue:.2f}")

    # Log dominant species and stop streaming if dominance threshold is met
    if dominant_percentage >= species_identification_percentage_dominance * 100:
        streaming_logger.info(f"Species {dominant_species} has been identified in {dominant_percentage:.2f}% of cases. Stopping streaming.")
        dbs = set(["ITS", "matK", "psbA-trnH", "rbcL"]) if db is None else ({db} if isinstance(db, str) else set(db))
        write_identification_results(output_dir, input_name, blastn_result_list[-1], dbs)
        return dominant_species, True
    return None, False

def run_streaming_basic_pipeline(input_name, input_fastq_path, output_dir, streaming_logger, db=None, windows=False, species_identification_percentage_dominance=80.0, block_size=500, minimum_block_amount_before_dominance_check=5, consensus_method="80_20_best_sequence", identification_method="blastn"):
    """
    Main function to run the streaming basic pipeline. It manages the overall workflow,
    logging, and tracking of time for different steps in the pipeline.
    """
    total_time_taken = 0
    total_time_taken_minimap2 = 0
    total_time_taken_racon = 0
    total_time_taken_blastn = 0

    total_reads = len(read_streaming_fastq(input_fastq_path))
    streaming_logger.info(f"Total number of reads: {total_reads}")

    blastn_result_list = []

    for i, block in enumerate(random_sampling_blocks(input_fastq_path, block_size, total_reads)):
        chunk_iteration_time_taken = 0

        streaming_logger.info(f"Running block {i+1}...")
        block_output_dir = setup_and_save_block(block, i+1, output_dir, streaming_logger)
        print(f"Block {i+1} saved at {block_output_dir}")
        block_fastq_path = os.path.join(block_output_dir, f"block_{i+1}.fastq")
        block_input_name = f"block_{i+1}"
 
        consensus_block_output_dir = os.path.join(block_output_dir, "consensus")
        os.makedirs(consensus_block_output_dir, exist_ok=True)

        consensus_time, minimap2_time, racon_time = run_consensus_for_block(block_input_name, block_fastq_path, consensus_method, block_output_dir, windows, streaming_logger)
        total_time_taken += consensus_time
        total_time_taken_minimap2 += minimap2_time
        total_time_taken_racon += racon_time
        chunk_iteration_time_taken += consensus_time

        blastn_result, blastn_time = run_identification_for_block(block_input_name, block_output_dir, db, identification_method, streaming_logger, consensus_block_output_dir, windows)
        total_time_taken += blastn_time
        total_time_taken_blastn += blastn_time
        chunk_iteration_time_taken += blastn_time
        blastn_result_list.append(blastn_result)

        # Log the result of the block for the identification
        for db, result in blastn_result.items():
            streaming_logger.info(f"Database: {db}")
            streaming_logger.info(f"Species: {result['species']}")
            streaming_logger.info(f"Alignment: {result['alignment']}")
            streaming_logger.info(f"E-value: {result['evalue']}")
            streaming_logger.info(f"--------------------------------------")

        streaming_logger.info(f"Block {i+1} completed. Total time taken: {chunk_iteration_time_taken:.2f} seconds.")
        streaming_logger.info(f"--------------------------------------")

        if i >= minimum_block_amount_before_dominance_check:
            dominant_species, goto_end = check_species_dominance(blastn_result_list, species_identification_percentage_dominance, streaming_logger, output_dir, input_name, db)
            if goto_end:
                break

    streaming_logger.info(f"Streaming pipeline completed. Total time taken: {total_time_taken:.2f} seconds.")
    streaming_logger.info(f"Total time taken for minimap2: {total_time_taken_minimap2:.2f} seconds.")
    streaming_logger.info(f"Total time taken for racon: {total_time_taken_racon:.2f} seconds.")
    streaming_logger.info(f"Total time taken for blastn: {total_time_taken_blastn:.2f} seconds.")
    total_time_taken_random_sampling = total_time_taken - total_time_taken_minimap2 - total_time_taken_racon - total_time_taken_blastn
    streaming_logger.info(f"Total time taken for random sampling: {total_time_taken_random_sampling:.2f} seconds.")
