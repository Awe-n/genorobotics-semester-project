from pipeline import pipeline_single
from lib.streaming.streaming_helpers.random_sampling import random_sampling_blocks, save_block_as_fastq
from lib.streaming.streaming_helpers.random_sampling import read_fastq

def run_streaming_basic_pipeline(input_name: str, input_fastq_path: str, output_dir: str, logger, wsl: bool = False):
    total_time_taken = 0
    total_time_taken_minimap2 = 0
    total_time_taken_racon = 0

    # # Run the pipeline
    # result = pipeline_single(input_fastq_filename, db, wsl)

    # # Print the result
    # print(result)

    # Calculate total number of reads
    total_reads = len(read_fastq(input_fastq_path))
    print(f"Total number of reads: {total_reads}")

    # Perform random sampling and save blocks
    block_size = 1000  # Set your block size
    for i, block in enumerate(random_sampling_blocks(input_fastq_path, block_size, total_reads)):
        # save_block_as_fastq(block, i+1, output_dir)
        print(f"Block {i+1} saved with {len(block)} reads")

    blocks = list(random_sampling_blocks(input_fastq_path, block_size, total_reads))

    # Sanity check
    if sanity_check(blocks):
        print("Sanity check passed: No duplicate reads found.")
    else:
        print("Sanity check failed: Duplicate reads detected.")


def sanity_check(blocks):
    """
    Check if there are any duplicate reads across the blocks.
    :param blocks: List of lists, where each inner list contains SeqRecord objects.
    :return: True if all reads are unique, False otherwise.
    """
    seen_ids = set()
    for block in blocks:
        for read in block:
            if read.id in seen_ids:
                return False  # Duplicate found
            seen_ids.add(read.id)
    return True

