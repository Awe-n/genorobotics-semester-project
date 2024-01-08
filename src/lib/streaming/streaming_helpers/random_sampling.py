import random
from Bio import SeqIO
import os

def read_fastq(fastq_filepath):
    with open(fastq_filepath, 'rt') as f:
        reads = list(SeqIO.parse(f, "fastq"))
    return reads

def random_sampling_blocks(fastq_file, block_size, total_reads):
    all_reads = read_fastq(fastq_file)
    all_indices = list(range(len(all_reads)))
    random.shuffle(all_indices)

    for i in range(0, total_reads, block_size):
        selected_indices = all_indices[i:i + block_size]
        yield [all_reads[idx] for idx in selected_indices]

def save_block_as_fastq(block, block_number, output_dir):
    output_file = os.path.join(output_dir, f"block_{block_number}.fastq")
    with open(output_file, "w") as output_handle:
        SeqIO.write(block, output_handle, "fastq")