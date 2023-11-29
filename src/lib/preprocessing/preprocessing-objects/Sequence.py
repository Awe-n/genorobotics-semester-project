import numpy as np
import warnings
from Bio import SeqRecord

from objects.PrimerAllAlignments import PrimerAllAlignments

class Sequence():

    def __init__(self, seq=None):
        self.sequence = None  # sequence [str]
        self.phred_scores = None  # Phred scores [np.ndarray]
        self.sequence_description = None  # sequence description from fastq files [str]
        self.mean_error_rate = None  # mean error rate computed from Phred scores [float]
        if seq is not None: self.set_sequence(seq)
        self.primer_forward_alignments = None # PrimerAllAlignments instance
        self.primer_reverse_alignments = None # PrimerAllAlignments instance
        self.primer_aligned = {"primer": None, "confidence": 0.0} # information about the primer aligned ("forward" or "reverse") and the primer alignment confidence associated
        self.primer_alignments_run = False # boolean specifying if run_primer_alignments() has been called
    
    def reset(self):
        """
        Reset the Sequence instance.
        :return: None
        """
        self.sequence = None
        self.phred_scores = None
        self.sequence_description = None
        self.mean_error_rate = None
        self.primer_forward_alignments = None
        self.primer_reverse_alignments = None
        self.primer_alignments_run = False
    
    def set_phred_scores(self, phred_scores):
        """
        Set Phred scores
        :param phred_scores: phred_scores [np.ndarray | 1D]
        :return: None
        """
        if type(self.phred_scores) != np.ndarray:
            phred_scores = np.array(phred_scores).reshape((-1,))
        self.phred_scores = phred_scores
        self.compute_mean_error_rate()
    
    def get_sequence_reversed(self):
        """
        Get the reverse sequence (5'-3' <-> 3'-5')
        /!\ Warning: reverse sequence has not the same meaning as the reverse primer
        :return: reversed sequence [str]
        """
        return self.sequence[::-1]
    
    def get_sequence_complementary(self):
        """
        Get the complementary sequence
        :return: complementary sequence [str]
        """
        return self.sequence.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()
    
    def get_sequence_variants(self):
        """
        Get the 4 variants of the sequence:
            - sequence (code 0|00)
            - complementary sequence (code 1|10)
            - reversed sequence (code 2|01)
            - reversed complementary sequence (code 3|11)
        :return: list of sequences [list of str | len=4]
        """
        sequence_comp = self.get_sequence_complementary()
        return [self.sequence, sequence_comp, self.sequence[::-1], sequence_comp[::-1]]
    
    def compute_mean_error_rate(self):
        """
        Compute and set the mean error rate (using the phred scores)
        :return: None
        """
        if self.phred_scores is not None:
            error_proba = 10 ** (-self.phred_scores / 10)
            self.mean_error_rate = error_proba.mean()
        else:
            self.mean_error_rate = None
    
    def flip_reverse(self):
        """
        Flip the Sequence to the reverse.
        :return: None
        """
        self.sequence = self.get_sequence_reversed()
        if self.phred_scores is not None:
            self.phred_scores = np.flip(self.phred_scores)
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.flip_reverse()
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.flip_reverse()
    
    def flip_complementary(self):
        """
        Flip the Sequence to the complementary.
        :return: None
        """
        self.sequence = self.get_sequence_complementary()
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.flip_complementary()
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.flip_complementary()
    
    def truncate(self, start, end=None):
        """
        Truncate the sequence.
        :param start: start of the truncation
        :param end: end of the truncation
        :return: None
        """
        if end is None: end = len(self.sequence)
        if start > len(self.sequence) or end > len(self.sequence):
            warnings.warn("Sequence.truncate(): invalid start/end given as arguments.")
        self.sequence = self.sequence[start:end]
        if self.phred_scores is not None:
            self.phred_scores = self.phred_scores[start:end]
        self.compute_mean_error_rate()
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.truncate(start, end)
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.truncate(start, end)
    
    def extend(self, start_added=0, end_added=0):
        """
        Extend the sequence with unknown nucleotides (N).
        :param start_added: number of nucleotides to add in the beginning [int]
        :param end_added: number of nucleotides to add at the end [int]
        :return: None
        """
        self.sequence = start_added*"N" + self.sequence + end_added*"N"
        if self.phred_scores is not None:
            self.phred_scores = np.concatenate((np.zeros(shape=(start_added,)), self.phred_scores, np.zeros(shape=(end_added,))))
        self.compute_mean_error_rate()
        if self.primer_forward_alignments is not None: self.primer_forward_alignments.extend(start_added, end_added)
        if self.primer_reverse_alignments is not None: self.primer_reverse_alignments.extend(start_added, end_added)
    
    def set_sequence(self, seq):
        """
        Set the sequence
        :param seq: sequence [Bio.SeqRecord.SeqRecord OR str]
        :return: None
        """
        if type(seq) is str:
            seq = ''.join([s for s in str(seq) if s in "ATGC"])
            if self.sequence is not None:
                found = False
                solvable = True
                for id, var in enumerate(self.get_sequence_variants()):
                    if seq in var:
                        if found: solvable = False # not possible if seq present in multiple var
                        found = True
                        loc_start = var.find(seq)
                        loc_end = loc_start + len(seq)
                        if var[loc_start+1:].find(seq) != -1: solvable = False # not possible if seq present multiple times in var
                        flip_complementary = id%2==1
                        flip_reverse = id//2==1
                        transformations_to_do = [flip_complementary, flip_reverse, loc_start, loc_end]
                if found == False or solvable == False:
                    self.reset()
                    self.sequence = seq
                else:
                    if transformations_to_do[0] == True: self.flip_complementary()
                    if transformations_to_do[1] == True: self.flip_reverse()
                    self.truncate(start=transformations_to_do[2], end=transformations_to_do[3])
            else:
                self.sequence = seq
        elif type(seq) is SeqRecord.SeqRecord:
            self.sequence = str(seq.seq)
            self.phred_scores = np.array(seq.letter_annotations['phred_quality'])
            self.sequence_description = str(seq.description)
        else:
            warnings.warn("The sequence type given to Sequence() class is not recognized.")
        self.compute_mean_error_rate()
    
    def set_sequence_with_best_primer_alignment(self):
        # Get the best directions from each primer alignments
        best_confidence_forward = self.primer_forward_alignments.get_best_direction()
        best_confidence_reverse = self.primer_reverse_alignments.get_best_direction()

        # Check if the best directions have been found
        error_forward, error_reverse = False, False
        if best_confidence_forward["confidence"] is None:
            best_confidence_forward["confidence"] = 0.0
            error_forward = True
        if best_confidence_reverse["confidence"] is None:
            best_confidence_reverse["confidence"] = 0.0
            error_reverse = True
        if error_forward and error_reverse:
            warnings.warn("[set_sequence_with_best_primer_alignments]: all primer alignments are empty")
            return None
        else:
            if best_confidence_forward["confidence"] > best_confidence_reverse["confidence"]:
                best_confidence = best_confidence_forward
                self.primer_aligned = {"primer": "forward", "confidence": best_confidence["confidence"]}
            else:
                best_confidence = best_confidence_reverse
                self.primer_aligned = {"primer": "reverse", "confidence": best_confidence["confidence"]}

        # Set the new sequence direction
        if best_confidence["complementary"]:
            self.flip_complementary()
        if best_confidence["reverse"]:
            self.flip_reverse()
    
    def truncate_sequence_using_best_primer_alignment(self):
        self.set_sequence_with_best_primer_alignment()
        if self.primer_aligned["primer"] == "forward":
            start = self.primer_forward_alignments.get_main_alignment_properties().get("start")
            self.truncate(start=start)
        elif self.primer_aligned["primer"] == "reverse":
            end = self.primer_reverse_alignments.get_main_alignment_properties().get("end")
            self.truncate(start=0, end=len(self.sequence)-end)
    
    def run_primer_alignments(self, primer_forward, primer_reverse):
        """
        Run primer_align for each Sequence.get_variants() for each primer (forward and reverse)
        This method initializes/sets the primers properties.
        :param primer_forward: 5'->3' forward primer [str]
                               ex: "ATGTCACCACAAACAGAGACTAAAGC"
        :param primer_reverse: 5'->3' reverse primer [str]
                               ex: "TCGCATGTACCTGCAGTAGC"
        :return: None
        """
        # Complementarize/reverse the primer_reverse (in order to put it in the same direction as primer_forward)
        primer_reverse = primer_reverse.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()
        primer_reverse = primer_reverse[::-1]

        self.primer_forward_alignments = PrimerAllAlignments(primer_forward)
        self.primer_reverse_alignments = PrimerAllAlignments(primer_reverse)
        self.primer_alignments_run = True
        for id, seq in enumerate(self.get_sequence_variants()):
            
            from helpers.alignment import align_primer
            
            ali_forward = align_primer(seq, primer_forward)
            ali_reverse = align_primer(seq, primer_reverse)
            reverse = id//2
            complementary = id%2
            self.primer_forward_alignments.set_alignment(ali_forward, complementary=complementary, reverse=reverse)
            self.primer_reverse_alignments.set_alignment(ali_reverse, complementary=complementary, reverse=reverse)
    
    def get_scores(self):
        self.compute_mean_error_rate()
        primer_forward_confidence = self.primer_forward_alignments.get_main_alignment_properties().get("confidence")
        primer_reverse_confidence = self.primer_reverse_alignments.get_main_alignment_properties().get("confidence")
        scores = {"mean_error_rate": self.mean_error_rate,
                  "primer_forward_confidence": primer_forward_confidence,
                  "primer_reverse_confidence": primer_reverse_confidence,
                  "primer_aligned": self.primer_aligned["primer"]}
        return scores
    
    def print_phred(self):
        """
        Print rescaled phred scores (between 0-9)
        This is useful to print it below a sequence
        -------------------
        score -> error rate |
            0 -> 25%-100%   |
            1 -> 6%-25%     |
            2 -> 1.5%-6%    |
            3 -> 0.4%-1.5%  |
            4 -> 0.1%-0.4%  |
        -------------------
        :return: None
        """
        # for i in range(10):
        #     print(10**((-6*i)/10))
        if self.phred_scores is None:
            print("phred_scores is None")
            return None
        phred_rescaled = np.floor(np.clip(self.phred_scores - 1, 0, 59) / 6).astype(int)  # phred scores rescaled between 0-9
        print("".join([str(phred_rescaled[i]) for i in range(phred_rescaled.size)]))
    
    def fastq_representation(self):
        """
        Return a fastq string representation of the Sequence instance
        :return: fastq string [str]
        """
        fastq_str = "@" + str(self.sequence_description) + "\n"
        fastq_str = fastq_str + self.sequence + "\n"
        fastq_str = fastq_str + "+\n"
        if self.phred_scores is not None:
            phred_scores = self.phred_scores
        else:
            phred_scores = len(self.sequence) * [0]
        fastq_str = fastq_str + "".join([chr(33 + c) for c in phred_scores]) + "\n"
        return fastq_str
    
    def __str__(self):
        """
        print(Sequence()) uses this function to format the Sequence() into string representation.
        Useful to print some informations about the Sequence() instance.
        :return: string representation of Sequence() instance [str]
        """
        string = "Sequence() instance\n"
        try:
            string = string + f"Mean error rate: {round(100 * self.mean_error_rate, 2)}%\n"
        except:
            pass
        try:
            primer_forward_alignment = self.primer_forward_alignments.get_formatted_main_alignment(self.sequence)
            properties_forward = self.primer_forward_alignments.get_main_alignment_properties()
            primer_forward_score, primer_forward_start, primer_forward_end = properties_forward.get("score"), properties_forward.get("start"), properties_forward.get("end")
            string = string + f"Primer forward: {primer_forward_alignment[2]} (score={round(primer_forward_score, 2)}, start={primer_forward_start}, end={primer_forward_end})\n"
            string = string + f"                {primer_forward_alignment[1]}\n"
            string = string + f"        {self.sequence[primer_forward_start - 8:primer_forward_start] + primer_forward_alignment[0] + self.sequence[primer_forward_end:primer_forward_end + 8]}\n"

            primer_reverse_alignment = self.primer_reverse_alignments.get_formatted_main_alignment(self.sequence)
            properties_reverse = self.primer_reverse_alignments.get_main_alignment_properties()
            primer_reverse_score, primer_reverse_start, primer_reverse_end = properties_reverse.get("score"), properties_reverse.get("start"), properties_reverse.get("end")
            string = string + f"Primer reverse: {primer_reverse_alignment[2]} (score={round(primer_reverse_score, 2)}, start={primer_reverse_start}, end={primer_reverse_end})\n"
            string = string + f"                {primer_reverse_alignment[1]}\n"
            string = string + f"        {self.sequence[primer_reverse_start - 8:primer_reverse_start] + primer_reverse_alignment[0] + self.sequence[primer_reverse_end:primer_reverse_end + 8]}\n"
        except:
            pass
        return string