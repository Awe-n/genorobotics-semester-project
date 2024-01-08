from pipeline import pipeline_single
from lib.streaming.streaming_helpers.random_sampling import random_sampling_blocks, save_block_as_fastq, read_streaming_fastq
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification
import os

from lib.general_helpers.configure_loggers import configure_consensus_logger
from lib.general_helpers.configure_loggers import configure_identification_logger

def run_streaming_basic_pipeline(input_name: str, input_fastq_path: str, output_dir: str, streaming_logger, db: str = None, wsl: bool = False):
    total_time_taken = 0
    total_time_taken_minimap2 = 0
    total_time_taken_racon = 0
    total_time_taken_blastn = 0

    goto_end = False

    # # Run the pipeline
    # result = pipeline_single(input_fastq_filename, db, wsl)

    # # Print the result
    # print(result)

    # Calculate total number of reads
    total_reads = len(read_streaming_fastq(input_fastq_path))
    streaming_logger.info(f"Total number of reads: {total_reads}")

    blastn_result_list = []

    # Perform random sampling and save blocks
    block_size = 500  # Set your block size: reads per block
    for i, block in enumerate(random_sampling_blocks(input_fastq_path, block_size, total_reads)):
        # Block Folder Setup
        block_output_dir = os.path.join(output_dir, f"block_{i+1}")
        os.makedirs(block_output_dir, exist_ok=True)

        # Save Block
        save_block_as_fastq(block, i+1, block_output_dir)
        streaming_logger.info(f"Block {i+1} saved with {len(block)} reads")
        block_fastq_path = os.path.join(block_output_dir, f"block_{i+1}.fastq")
        block_input_name = f"block_{i+1}"

        # Consensus Block
        consensus_block_output_dir = os.path.join(block_output_dir, "consensus")
        os.makedirs(consensus_block_output_dir, exist_ok=True)

        consensus_logger = configure_consensus_logger(consensus_block_output_dir, block_input_name)
        consensus_logger.info(f"Logging set up at {consensus_block_output_dir}/{block_input_name}_consensus_pipeline_log.log")
        streaming_logger.info(f"Running consensus pipeline for block {i+1}...")

        consensus_method = "80_20_best_sequence"
        consensus_total_time_taken, consensus_total_time_taken_minimap2, consensus_total_time_taken_racon = run_consensus(block_input_name, block_fastq_path, consensus_method, consensus_block_output_dir, consensus_logger, wsl)
        total_time_taken += consensus_total_time_taken
        total_time_taken_minimap2 += consensus_total_time_taken_minimap2
        total_time_taken_racon += consensus_total_time_taken_racon

        streaming_logger.info(f"Consensus pipeline for block {i+1} completed")

        # Identification Block
        identification_block_output_dir = os.path.join(block_output_dir, "blastn")
        os.makedirs(identification_block_output_dir, exist_ok=True)

        # /block_1/block_1.fastq/block_1_final_consensus.fasta'

        identification_logger = configure_identification_logger(identification_block_output_dir, block_input_name)
        identification_logger.info(f"Logging set up at {identification_block_output_dir}/{block_input_name}_identification_pipeline_log.log")
        streaming_logger.info(f"Running identification pipeline for block {i+1}...")

        blastn_result, identification_total_time_taken_blastn = run_identification(input_name=block_input_name, expedition_name=None, input_path=consensus_block_output_dir, output_dir=identification_block_output_dir, db=db, logger=identification_logger)
        total_time_taken += identification_total_time_taken_blastn
        total_time_taken_blastn += identification_total_time_taken_blastn

        streaming_logger.info(f"Identification pipeline for block {i+1} completed")
        streaming_logger.info(f"Result for block {i+1}:")
        
        for db, result in blastn_result.items():
            streaming_logger.info(f"Database: {db}")
            streaming_logger.info(f"Species: {result['species']}")
            streaming_logger.info(f"Alignment: {result['alignment']}")
            streaming_logger.info(f"E-value: {result['evalue']}")
            streaming_logger.info(f"--------------------------------------")

        blastn_result_list.append(blastn_result)

        # if since the beginning of the streaming, the same species has been identified 5 times (not necessarily in a row), then stop the streaming and log the species
        if i >= 4:
            species_list = []
            for blastn_result in blastn_result_list:
                for db, result in blastn_result.items():
                    species_list.append(result['species'])
            species_count = {}
            for species in species_list:
                if species in species_count:
                    species_count[species] += 1
                else:
                    species_count[species] = 1
            for species, count in species_count.items():
                if count >= 5:
                    streaming_logger.info(f"Species {species} identified 5 times. Stopping the streaming.")
                    goto_end = True
                
        if goto_end:
            break
    
    streaming_logger.info(f"Streaming pipeline completed. Total time taken: {total_time_taken:.2f} seconds.")
    streaming_logger.info(f"Total time taken for minimap2: {total_time_taken_minimap2:.2f} seconds.")
    streaming_logger.info(f"Total time taken for racon: {total_time_taken_racon:.2f} seconds.")
    streaming_logger.info(f"Total time taken for blastn: {total_time_taken_blastn:.2f} seconds.")
    total_time_taken_random_sampling = total_time_taken - total_time_taken_minimap2 - total_time_taken_racon - total_time_taken_blastn
    streaming_logger.info(f"Total time taken for random sampling: {total_time_taken_random_sampling:.2f} seconds.")
