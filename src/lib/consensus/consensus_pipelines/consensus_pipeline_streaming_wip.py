import os
from Bio import SeqIO
from lib.general_helpers.run_command import run_command
import logging
from Bio.Seq import Seq

def invert_bases(seq):
    """Inverts the DNA bases in a sequence."""
    valid_bases = set('ATCGatcg')
    if all(base in valid_bases for base in seq):
        return seq.translate(str.maketrans('ATCGatcg', 'TAGCtagc'))
    else:
        # Handle sequences with invalid characters (like 'N')
        logging.warning(f"Sequence contains invalid characters and will not be inverted: {seq}")
        return seq

def refine_consensus_sequence(consensus_path, output_dir, block_num):
    updated_sequences = []
    try:
        sequences = list(SeqIO.parse(consensus_path, 'fasta'))
        first_bases = [str(seq.seq)[0].upper() for seq in sequences if len(seq) > 0]

        at_count = first_bases.count('A') + first_bases.count('T')
        cg_count = first_bases.count('C') + first_bases.count('G')
        total = at_count + cg_count

        most_common = 'AT' if at_count >= cg_count else 'CG'
        least_common = 'CG' if most_common == 'AT' else 'AT'
        percentage_dominance = (at_count / total * 100) if most_common == 'AT' else (cg_count / total * 100)

        logging.info(f"Percentage Dominance of {most_common} at First Position: {percentage_dominance:.2f}%")

        for seq in sequences:
            if seq.seq[0] in least_common:
                inverted_seq_str = invert_bases(str(seq.seq))
                seq.seq = Seq(inverted_seq_str)
            updated_sequences.append(seq)

        updated_consensus_path = os.path.join(output_dir, f"consensus_updated_{block_num}.fasta")
        with open(updated_consensus_path, 'w') as output_file:
            SeqIO.write(updated_sequences, output_file, 'fasta')

        return updated_consensus_path
    except Exception as e:
        logging.error(f"Error processing file {consensus_path}: {e}")
        return None
    


def process_block(reads, output_dir, block_num, previous_consensus=None):
    block_fastq_path = os.path.join(output_dir, f"block_{block_num}.fastq")
    SeqIO.write(reads, block_fastq_path, "fastq")

    alignment_path = os.path.join(output_dir, f"alignment_{block_num}.paf")
    consensus_path = os.path.join(output_dir, f"consensus_{block_num}.fasta")

    # Align reads (to each other or to previous consensus)
    if previous_consensus:
        alignment_command = f"minimap2 -x map-ont {previous_consensus} {block_fastq_path} > {alignment_path}"
    else:
        alignment_command = f"minimap2 -x ava-ont {block_fastq_path} {block_fastq_path} > {alignment_path}"
    _, _ = run_command(alignment_command)

    # Generate consensus sequence with racon
    racon_command = f"racon -m 8 -x -6 -g -8 -w 500 {block_fastq_path} {alignment_path} {block_fastq_path} > {consensus_path}"
    _, _ = run_command(racon_command)

    # Refine the consensus sequence
    updated_consensus_path = refine_consensus_sequence(consensus_path, output_dir, block_num)

    return updated_consensus_path



def run_consensus_streaming_pipeline_wip(input_name: str, input_fastq_path: str, output_dir: str, block_size=3000):
    os.makedirs(output_dir, exist_ok=True)
    all_consensus_sequences = []
    block_number = 0
    current_block_reads = []

    logging.info("Processing blocks...")

    with open(input_fastq_path, 'r') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            current_block_reads.append(record)
            
            if len(current_block_reads) == block_size:
                # Process the current block to get a consensus sequence
                previous_consensus = all_consensus_sequences[-1] if all_consensus_sequences else None
                consensus_path = process_block(current_block_reads, output_dir, block_number, previous_consensus)
                all_consensus_sequences.append(consensus_path)

                # Reset for the next block
                current_block_reads = []
                block_number += 1

                # if block_number == 2 then stop the function and return the consensus_path
                if block_number == 2:
                    return consensus_path

        # Process the last block if it's not empty
        if current_block_reads:
            previous_consensus = all_consensus_sequences[-1] if all_consensus_sequences else None
            consensus_path = process_block(current_block_reads, output_dir, block_number, previous_consensus)
            all_consensus_sequences.append(consensus_path)

    # TODO: Implement ranking and selection of top consensus sequences
    # This might involve parsing the consensus sequences, evaluating their quality, length, etc.

    return all_consensus_sequences[:4]  # Return the top 4 consensus sequences