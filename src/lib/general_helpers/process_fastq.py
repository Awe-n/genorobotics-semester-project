import glob, gzip, math, random
import numpy as np
from Bio import SeqIO, Align
import os
import gzip
import shutil



## AWEN'S CODE

def split_fastq(input_path, output_dir, base_name, percentile=20):
    """
    Split a FASTQ file into two separate files based on sequence length.

    Args:
        input_path (str): Path to the input FASTQ file.
        output_dir (str): Directory where the output files will be saved.
        base_name (str): Base name for the output files.
        percentile (int, optional): The percentile value to split the sequences. Defaults to 20.

    Returns:
        tuple: A tuple containing the paths to the top percentile sequences file and the remaining sequences file.
    """
    # Read sequences and sort by length
    sequences = list(SeqIO.parse(input_path, "fastq"))
    sequences.sort(key=lambda x: len(x), reverse=True)

    # Split into top percentile and remaining sequences
    split_index = len(sequences) * percentile // 100
    top_sequences = sequences[:split_index]
    remaining_sequences = sequences[split_index:]

    # Write to separate files
    top_sequences_path = os.path.join(output_dir, f"{base_name}_top{percentile}.fastq")
    remaining_sequences_path = os.path.join(output_dir, f"{base_name}_remaining{100-percentile}.fastq")
    SeqIO.write(top_sequences, top_sequences_path, "fastq")
    SeqIO.write(remaining_sequences, remaining_sequences_path, "fastq")

    return top_sequences_path, remaining_sequences_path


## MILOU'S CODE

def extract_gz(src_file, dst_file):
    """
    Extracts a gzipped file.

    Args:
        src_file (str): Path to the source gzipped file.
        dst_file (str): Path to the destination file.

    Returns:
        None
    """
    with gzip.open(src_file, 'rb') as f_in:
        with open(dst_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)






## JEREMY'S CODE

def read_fastq(fastq_filepath=None):
    """
    Read a fastq file and return the reads.
    :param fastq_filepath: filepath of the .fastq file [str]
    :return: reads from the fastq file [list of Bio.SeqRecord.SeqRecord]
    """
    if fastq_filepath is None: fastq_filepath = "data/rbcL_Qiagen_tomato.fastq" # default path (example)
    if fastq_filepath.lower().endswith('.gz'):
        f = gzip.open(fastq_filepath, 'rt')
    else:
        f = open(fastq_filepath, 'rt')
    reads = []
    for read in SeqIO.parse(f, "fastq"):
        reads.append(read)
    return reads

def concatenate_fastq(src_folder=None, dst=None):
    """
    Concatenate all .fastq from a folder into a single .fastq file.
    :param folder: folder containing the .fastq files (usually fastq_pass folder) [str]
    :param dst: destination file (.fastq) [str]
    :return: None
    """
    if src_folder is None: src_folder = "fastq_pass" # default folder (example)
    if dst is None: dst = f"{src_folder}/concatenation.fastq" # default destination (example)
    if dst[-6:] != ".fastq": dst = f"{dst}.fastq"
    def get_file_iterator(filename):
        if filename.lower().endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'rt')
        return f
    fastq_list = [fastq for fastq in glob.glob(f"{src_folder}/*.fastq*")]
    fastq_iterators = [SeqIO.parse(get_file_iterator(fastq), "fastq") for fastq in fastq_list]
    while True:
        for fq in fastq_iterators:
            try:
                SeqIO.write(next(fq), open(dst,"at"), "fastq")
            except StopIteration:
                fastq_iterators.remove(fq)
        if len(fastq_iterators) == 0:
            break

def get_substring(string, k): # ACTUALLY NOT USED IN THE CODE
    """
    Get substring of length k from a string.
    :param string: original string [str]
    :param k: length of substrings [int]
    :return: list of substrings [list of str | len=len(string)-k+1]
    """
    return [string[i:i+k] for i in range(0, len(string)-k+1)]

def get_subprimer(primer, k, error_ratio=0): # ACTUALLY NOT USED IN THE CODE
    """
    Get subprimers (substrings of a primer) of length k and the
    :param primer: primer sequence [str]
    :param k: length of subprimers [int]
    :param error_ratio: error_ratio allowed [float]
    :return: list of subprimers (subprimer sequence, left characters to try, right characters to try) [list of [str, int, int]]
    """
    ratio = 1.0 + error_ratio
    return [[primer[i:i+k], math.ceil(ratio*i), math.ceil(ratio*(len(primer)-i-k))] for i in range(0, len(primer)-k+1)]

def create_random_sequence(length=500, seed=None):
    """
    Generate a random sequence.
    :param length: length of the sequence generated [int]
    :param seed: random seed [int]
    :return: random sequence [str]
    """
    if seed is not None: random.seed(seed)
    return "".join(["ATGC"[random.randint(0,3)] for i in range(length)])

def damage_sequence(sequence, mutation_rate=0.05, deletion_rate=0.05, insertion_rate=0.05):
    """
    Damage a sequence randomly with mutation (substitution), deletion and insertion
    Warning: to avoid infinite loop, don't set the insertion_rate to 1.0 (or near)
    :param sequence: original sequence to damage [str]
    :param mutation_rate: mutation rate (between 0.0 and 1.0) [float]
    :param deletion_rate: deletion_rate (between 0.0 and 1.0) [float]
    :param insertion_rate: insertion_rate (between 0.0 and 1.0) [float]
    :return: damaged sequence [str]
    """
    if mutation_rate < 0 or mutation_rate > 1:
        raise Exception("[damage_sequence] mutation_rate is incorrect (must be between 0.0 and 1.0)")
    if deletion_rate < 0 or deletion_rate > 1:
        raise Exception("[damage_sequence] deletion_rate is incorrect (must be between 0.0 and 1.0)")
    if insertion_rate < 0 or insertion_rate > 1:
        raise Exception("[damage_sequence] insertion_rate is incorrect (must be between 0.0 and 1.0)")
    sequence = "".join(["ATGC".replace(b, '')[random.randint(0, 2)] if random.random() < mutation_rate else b for b in sequence]) # mutation / substitution
    sequence = "".join(['' if random.random() < deletion_rate else b for b in sequence]) # deletion
    insertion_extension_rate = insertion_rate  # can be changed if the extension rate is different
    def get_insert(extension_rate=insertion_extension_rate):
        insert = "ATGC"[random.randint(0,3)]
        while random.random() < extension_rate: insert += "ATGC"[random.randint(0,3)] # extension (after insertion)
        return insert
    sequence = "".join([sequence[i:i+1] + get_insert() if random.random() < insertion_rate else sequence[i:i+1] for i in range(len(sequence)+1)]) # insertion
    return sequence

def sigmoid_primer(x):
    # temporary function (ASUP later) -> used to give a probability to a primer score
    x = 14 * x - 7
    if x >= 100:
        sig = 1.0
    elif x <= -100:
        sig = 0.0
    else:
        sig = 1/(1+np.exp(-x))
    return sig