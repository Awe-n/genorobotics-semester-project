import random
from Bio import SeqIO
import os

def read_streaming_fastq(fastq_filepath):
    with open(fastq_filepath, 'rt') as f:
        reads = list(SeqIO.parse(f, "fastq"))
    return reads

def random_sampling_blocks(fastq_file, block_size, total_reads):
    all_reads = read_streaming_fastq(fastq_file)
    all_indices = list(range(len(all_reads)))
    random.shuffle(all_indices)

    for i in range(0, total_reads, block_size):
        selected_indices = all_indices[i:i + block_size]
        yield [all_reads[idx] for idx in selected_indices]

def save_block_as_fastq(block, block_number, output_dir):
    output_file = os.path.join(output_dir, f"block_{block_number}.fastq")
    with open(output_file, "w") as output_handle:
        SeqIO.write(block, output_handle, "fastq")

def random_sampling_sanity_check(blocks):
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