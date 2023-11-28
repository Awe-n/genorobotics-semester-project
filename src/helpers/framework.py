"""
Tester les algorithmes pour le choix d'un algorithme efficace et juste (vitesse et justesse)

Pour augmenter la rapidite, utiliser la parallelisation?

For now, we are doing alignments after concatenating the reads we've got from NP

input: 
    - sequence retrieved from a chosen algorithm
    - Reference sequence
output:
    - alignment score

"""

# from sequence_preprocessing import *
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from main.sequence_preprocessing import *
from itertools import chain

from helpers.fastq import read_fastq

def load_fasta(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq = seq_record.seq
    return seq

class ConsensusSeq:
    def __init__(self, name, filename=None, seq=None):
        self.name = name  # name of the algorithm
        self.filename = None
        if filename:
            self.filename = filename
            self.sequence = read_fastq(filename)
            # self.sequence = load_fasta(filename)  # output in fastq?
        elif seq:
            self.sequence = seq
        self.length = len(self.sequence)
        # metrics
        self.score = 0  # alignment score
        self.evalue = 0
        self.coverage = 0  # percentage of coverage

        # Example of consensus? https://stackoverflow.com/questions/73702044/how-to-get-a-consensus-of-multiple-sequence-alignments-using-biopython

    def get_consensus_seq(self):
        common_alignment = MultipleSeqAlignment(
            chain(*AlignIO.parse(self.filename, "fastq"))
        )
        summary = SummaryInfo(common_alignment)
        consensus = summary.dumb_consensus(0.7, "N")
        self.sequence = consensus

    def align_ref(self, ref):
        aligner = Bio.Align.PairwiseAligner()
        alignments = aligner.align(self.sequence, ref)
        self.score = alignments.score
        self.coverage = self.length / len(ref)
        print("Score: " + str(self.score))
        print("Coverage: " + str(self.coverage))


ref = load_fasta("rbcL_NCBI_tomato.fasta")
output = ConsensusSeq("test", "BestSequences_rbcL_Qiagen_tomato_5000.fastq")
# print(MultiSeq.consensus)
output.get_consensus_seq()
output.align_ref(ref)
