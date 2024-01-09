import random
from Bio import SeqIO
import os

def read_streaming_fastq(fastq_filepath):
    """
    Reads a FASTQ file and returns a list of SeqRecord objects representing the reads.
    :param fastq_filepath: Path to the FASTQ file.
    :return: A list of SeqRecord objects.
    """
    with open(fastq_filepath, 'rt') as f:
        reads = list(SeqIO.parse(f, "fastq"))
    return reads

def random_sampling_blocks(fastq_file, block_size, total_reads):
    """
    Generates blocks of reads of a specified size from the FASTQ file in a random manner.
    :param fastq_file: Path to the FASTQ file.
    :param block_size: The number of reads per block.
    :param total_reads: Total number of reads to consider for generating blocks.
    :yield: A block of reads (list of SeqRecord objects).
    """
    all_reads = read_streaming_fastq(fastq_file)
    all_indices = list(range(len(all_reads)))
    random.shuffle(all_indices)

    for i in range(0, total_reads, block_size):
        selected_indices = all_indices[i:i + block_size]
        yield [all_reads[idx] for idx in selected_indices]

def save_block_as_fastq(block, block_number, output_dir):
    """
    Saves a block of reads as a FASTQ file.
    :param block: A list of SeqRecord objects to be saved.
    :param block_number: The block number (used in naming the output file).
    :param output_dir: Directory where the FASTQ file will be saved.
    """
    output_file = os.path.join(output_dir, f"block_{block_number}.fastq")
    with open(output_file, "w") as output_handle:
        SeqIO.write(block, output_handle, "fastq")

def random_sampling_sanity_check(blocks):
    """
    Performs a sanity check to ensure no duplicate reads exist across blocks.
    :param blocks: A list of blocks, each containing SeqRecord objects.
    :return: Boolean value indicating whether all reads are unique (True) or if duplicates exist (False).
    """
    seen_ids = set()
    for block in blocks:
        for read in block:
            if read.id in seen_ids:
                return False  # Duplicate found
            seen_ids.add(read.id)
    return True