import warnings

class PrimerAllAlignments():

    def __init__(self, primer):
        self.primer = primer  # primer sequence [str]
        self.alignments = [[None, None], [None, None]] # 4 instances of PrimerAlignment [[normal, reverse], [complementary, reverse_complementary]]

    def set_alignment(self, PrimerAlign, complementary=False, reverse=False):
        """
        Setter of a PrimerAlignment
        :param PrimerAlign: PrimerAlignment to set [PrimerAlignment instance]
        :param complementary: is the alignment for complementary sequence [bool]
        :param reversed: is the alignment for reversed sequence [bool]
        :return: None
        """
        PrimerAlign.sequence_complementary = complementary
        PrimerAlign.sequence_reverse = reverse
        self.alignments[int(complementary)][int(reverse)] = PrimerAlign
        
    def get_formatted_main_alignment(self, sequence):
        if self.alignments[0][0] is not None:
            main_primer_alignment = self.alignments[0][0].get_formatted_alignment(sequence, self.primer)
        else:
            main_primer_alignment = None
            warnings.warn("PrimerAllAlignments.get_formatted_main_alignment(): the main primer_alignment has not been set. It can't be accessed.")
        return main_primer_alignment
    
    def get_main_alignment_properties(self):
        if self.alignments[0][0] is not None:
            properties = {"score": self.alignments[0][0].primer_score,
                          "start": self.alignments[0][0].primer_start,
                          "end": self.alignments[0][0].primer_end,
                          "confidence": self.alignments[0][0].primer_confidence}
            return properties
        else:
            warnings.warn("PrimerAllAlignments.get_main_alignment_properties(): the main primer_alignment has not been set. It can't be accessed.")
            return dict()
        
    def flip_reverse(self):
        new_alignments = [[self.alignments[0][1], self.alignments[0][0]], [self.alignments[1][1], self.alignments[1][0]]]
        self.alignments = new_alignments
        for id1 in range(2):
            for id2 in range(2):
                self.alignments[id1][id2].sequence_reverse = not self.alignments[id1][id2].sequence_reverse

    def flip_complementary(self):
        new_alignments = [[self.alignments[1][0], self.alignments[1][1]], [self.alignments[0][0], self.alignments[0][1]]]
        self.alignments = new_alignments
        for id1 in range(2):
            for id2 in range(2):
                self.alignments[id1][id2].sequence_complementary = not self.alignments[id1][id2].sequence_complementary

    def truncate(self, start, end):
        for id1 in range(2):
            for id2 in range(2):
                if self.alignments[id1][id2].sequence_reverse:
                    start_local, end_local = end, start
                else:
                    start_local, end_local = start, end
                if start_local > self.alignments[id1][id2].primer_start:
                    self.alignments[id1][id2] = None
                else:
                    self.alignments[id1][id2].primer_start = self.alignments[id1][id2].primer_start - start_local
                    self.alignments[id1][id2].primer_end = self.alignments[id1][id2].primer_end - start_local

    def extend(self, start_added, end_added):
        for id1 in range(2):
            for id2 in range(2):
                if self.alignments[id1][id2].sequence_reverse:
                    start_added_local, end_added_local = end_added, start_added
                else:
                    start_added_local, end_added_local = start_added, end_added
                self.alignments[id1][id2].primer_start = self.alignments[id1][id2].primer_start + start_added_local
                self.alignments[id1][id2].primer_end = self.alignments[id1][id2].primer_end + start_added_local

    def get_best_direction(self):
        best = {"confidence": None, "complementary": None, "reverse": None}
        for is_complementary in range(2):
            for is_reverse in range(2):
                primer_confidence = self.alignments[is_complementary][is_reverse].primer_confidence
                if primer_confidence is not None:
                    if best["confidence"] is None:
                        best["confidence"] = 0.0
                    if primer_confidence > best["confidence"]:
                        best["confidence"] = primer_confidence
                        best["complementary"] = bool(is_complementary)
                        best["reverse"] = bool(is_reverse)
        return best