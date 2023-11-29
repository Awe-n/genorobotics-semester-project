from helpers.alignment import align_primer

class PrimerAlignment():

    def __init__(self, primer_score=None, primer_start=None, primer_end=None, primer_confidence=None, primer_alignment=None, sequence_complementary=False, sequence_reverse=False):
        self.primer_score = primer_score  # alignment score between primer and sequence [float]
        self.primer_start = primer_start  # start position of primer in sequence [int]
        self.primer_end = primer_end  # end position of primer in sequence [int]
        self.primer_confidence = primer_confidence
        self.primer_alignment = primer_alignment  # alignment information formatted (between primer and sequence) [list of str]
        self.sequence_complementary = sequence_complementary
        self.sequence_reverse = sequence_reverse

    def get_formatted_alignment(self, sequence, primer):
        """
        Get formatted strings of the primer alignment
        :param sequence: sequence [str]
                         the sequence given must be "reverse" if self.sequence_reverse = True
                         the sequence given must be "complementary" if self.sequence_complementary = True
        :return: list of 3 strings containing formatted sequence alignment (truncated sequence, alignment, primer) [list of 3 str]
        """
        if self.sequence_reverse: sequence = sequence[::-1]
        if self.sequence_complementary: sequence = sequence.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()
        seq_primed = sequence[self.primer_start:self.primer_end]
        return align_primer(seq_primed, primer, return_formatted_alignment=True).primer_alignment
    
