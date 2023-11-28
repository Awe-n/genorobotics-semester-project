from helpers.fastq import sigmoid_primer, create_random_sequence
import numpy as np
import random

from Bio import Align

def align_primer(sequence, primer, aligner_scores=None, print_alignment=False, return_formatted_alignment=False):
    """
    Global alignment of primer inside sequence.
    :param sequence: sequence where the primer is searched [str]
    :param primer: primer searched [str]
    :param print_alignment: print alignment results [bool]
    :return: primer alignment results [float, int, int, [str, str, str]]
            score: alignment score [float]
             loc_start: start position of the primer in the sequence [int]
             loc_end: end position of the primer in the sequence [int]
             primer_alignment: list of 3 strings containing formatted sequence alignment (truncated sequence, alignment, primer) [list of 3 str]
    """
    # Aligner settings
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    # print(Align.substitution_matrices.load()) # print possible substitution matrices
    if aligner_scores is None: aligner_scores = dict()
    aligner.match_score = aligner_scores.get("match_score") if aligner_scores.get("match_score") is not None else 1 #3
    aligner.mismatch_score = aligner_scores.get("mismatch_score") if aligner_scores.get("mismatch_score") is not None else -1 #-2
    aligner.gap_score = aligner_scores.get("gap_score") if aligner_scores.get("gap_score") is not None else -2
    aligner.query_end_gap_score = 0 # must be 0 (to allow gap around the primer for global alignment) [query ~ sequence here]
    aligner.target_end_gap_score = -100 # set to -100 if the primer should be fully integrated in the sequence [target ~ primer here]

    # Alignment
    alignments = aligner.align(sequence, primer)
    # Get the best alignment only
    alignment = max(alignments, key=lambda x: x.score) # old version (slower): # alignment = sorted(alignments)[0]
    if print_alignment:
        print(f"Score = {round(alignment.score,1)}")
        print(alignment)

    # 'score', 'loc_start', 'loc_end' computations
    score = alignment.score/len(primer)
    loc_start = alignment.aligned[0][0][0]
    loc_end = alignment.aligned[0][-1][1]

    # 'primer_alignment' contains formatted string used to print the alignment results
    primer_alignment = None
    if return_formatted_alignment:
        primer_alignment = []
        aform = alignment.format() # format() calls _format_pretty()
        aform_linebreak1 = aform.find('\n')
        aform_linebreak2 = aform.find('\n', aform_linebreak1 + 1)
        aform_linebreak3 = aform.find('\n', aform_linebreak2 + 1)
        primer_position_str = aform[aform_linebreak2+1:aform_linebreak3].replace("A", "N").replace("T", "N").replace("G", "N").replace("C", "N")
        pos1 = primer_position_str.find("N")
        pos2 = primer_position_str.rfind("N") + 1
        primer_alignment.append(aform[0:aform_linebreak1][pos1:pos2])
        primer_alignment.append(aform[aform_linebreak1+1:aform_linebreak2][pos1:pos2])
        primer_alignment.append(aform[aform_linebreak2+1:aform_linebreak3][pos1:pos2])

    conf = sigmoid_primer(score)
    
    from objects.PrimerAlignment import PrimerAlignment

    return PrimerAlignment(primer_score=score, primer_start=loc_start, primer_end=loc_end, primer_confidence=conf, primer_alignment=primer_alignment)

def align_primer_random_scores(seq_length=500, primer_length=26, iterations=1000, seed=0, print_scores=True):
    """
    Evaluate alignment scores of random alignments (useful to determine score threshold to avoid the use of random alignments)
    :param seq_length: length of sequence [int]
    :param primer_length: length of primer [int]
    :param iterations: number of iterations (more iterations -> more statistically significant but slower) [int]
    :param seed: random seed [int]
    :param print_scores: print the results/scores [bool]
    :return: 99.9% percentile score [float]
    """
    if seed is not None: random.seed(seed)
    scores = np.empty(shape=(iterations,))
    for it in range(iterations):
        primer = create_random_sequence(primer_length)
        seq = create_random_sequence(seq_length)
        scores[it] = align_primer(seq, primer).primer_score
        # if it%1000==0: print(it)
    scores = scores/np.mean(scores)
    if print_scores:
        print(f"Random alignments scores (seq_length={seq_length}, primer_length={primer_length}, it={iterations}):")
        print(f"\tMean: {np.mean(scores)}")
        print(f"\tStd: {np.std(scores)}")
        for percent in [50, 75, 90, 95, 99, 99.9]:
            print(f"\tBest {round(100-percent,1)}%: {round(np.percentile(scores, percent),2)}")
    return scores