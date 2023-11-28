import os
from Bio import SeqIO
from Bio import pairwise2

def compute_alignment_score(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq, one_alignment_only=True)
    if alignments:
        alignment = alignments[0]
        return alignment.score
    else:
        return 0  # No alignment found, score is 0
    
def select_best_alignment(sequences, query_sequence):
    best_sequence = None
    best_score = 0

    for sequence in sequences:
        alignment_score = compute_alignment_score(sequence, query_sequence)  # Implement this function to compute alignment score
        if alignment_score > best_score:
            best_score = alignment_score
            best_sequence = sequence

    return best_sequence