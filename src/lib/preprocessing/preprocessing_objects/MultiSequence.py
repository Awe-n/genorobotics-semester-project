import warnings, time, copy, pickle
import numpy as np
from helpers.fastq import read_fastq
from objects.Sequence import Sequence

class MultiSequence():
    def __init__(self, filepath=None):
        self.sequences = [] # list of Sequence() instances [list of Sequence()]
        self.filepath = filepath # filepath used to load data [str]
        self.consensus = None # consensus sequence [str]
        self.consensus_phred_scores = None # consensus sequence Phred scores [np.ndarray]
        self.run_infos = {"primer_alignments_run": False, # set primer_alignments_run to True when run_primer_alignments() has been run [bool]
                          "ranking_computed": False}
        self.ranking_forward = None
        self.ranking_reverse = None
        if self.filepath is not None:
            self.read_file()
    
    def read_file(self):
        """
        Read the file in self.filepath, 2 file extension are accepted:
            - .fast: load the data
            - .pickle: load the MultiSequence() previously saved
        /!\ Warning: This function is intended to be used once (for a MultiSequence() instance)
                     It means that multiple .fastq files should be concatenated before being read by this method
        :return: None
        """
        if self.filepath is not None:
            extension = self.filepath[self.filepath.rfind('.'):]
            if extension == ".fastq":
                reads = read_fastq(fastq_filepath=self.filepath)
                for read in reads:
                    self.sequences.append(Sequence(read))
            elif extension == ".pickle":
                self.load(self.filepath)
            else:
                warnings.warn("Error in MultiSequence(): filepath extension not recognised during read_file().")
        else:
            warnings.warn("Error in MultiSequence(): no filepath found during read_file().")
    
    def write_to_fastq(self, dst=None):
        """
        Write all the sequences in a single .fastq formatted file
        :param dst: filepath of the fastq file destination [str]
        :return: None
        """
        if dst is None: dst = "MultiSequence_saved.fastq"  # default destination (example)
        if dst[-6:] != ".fastq": dst = f"{dst}.fastq"
        with open(dst, mode="w") as f:
            for Seq in self.sequences:
                f.write(Seq.fastq_representation())
    
    def save(self, dst=None):
        """
        Save the current instance of MultiSequence() to the file dst
        :param dst: filepath of destination [str]
        :return: None
        """
        if dst is None: dst = "MultiSequence_saved.pickle"  # default destination (example)
        if dst[-7:] != ".pickle": dst = f"{dst}.pickle"
        pickle.dump(self.__dict__, open(dst, mode='wb'))
    
    def load(self, src=None):
        """
        Load a MultiSequence() instance saved
        :param src: filepath of the pickle file source [str]
        :return: None
        """
        if src is None: src = "MultiSequence_saved.pickle"  # default source (example)
        if src[-7:] != ".pickle": src = f"{src}.pickle"
        state = pickle.load(open(src, mode='rb'))
        self.__dict__.update(state)
    
    def run_primer_alignments(self, primer_forward, primer_reverse):
        """
        Run primer_align on each self.sequences for each Sequence.get_variants() for each primer (forward and reverse)
        This method initializes/sets the primers properties in each Sequence()
        This method prints a progress bar.
        :param primer_forward: 5'->3' forward primer [str]
                               ex: "ATGTCACCACAAACAGAGACTAAAGC"
        :param primer_reverse: 5'->3' reverse primer [str]
                               ex: "TCGCATGTACCTGCAGTAGC"
        :return: None
        """
        time_last = time.time()
        mean_t = 0.0
        for Seq_id in range(len(self.sequences)):
            self.sequences[Seq_id].run_primer_alignments(primer_forward, primer_reverse)
            self.sequences[Seq_id].set_sequence_with_best_primer_alignment()

            delta_t = time.time() - time_last
            time_last = time.time()
            if mean_t == 0: mean_t = delta_t + 1e-10
            mean_t = mean_t * 0.99 + delta_t * 0.01 # "running" time
            print(f"\rrun_primer_alignments progress: {Seq_id+1}/{len(self.sequences)} ({round(1 / mean_t)} reads/sec)", end='')
        print("")
        self.run_infos["primer_alignments_run"] = True
   
    def compute_sequence_ranking(self):
        sequences_scores = []
        for Seq_id in range(len(self.sequences)):
            scores = self.sequences[Seq_id].get_scores()
            sequences_scores.append((Seq_id, scores))
        ranking_forward = sorted(sequences_scores, key=lambda x: (x[1]["primer_forward_confidence"], x[1]["mean_error_rate"]), reverse=True)
        self.ranking_forward = [r[0] for r in ranking_forward]
        ranking_reverse = sorted(sequences_scores, key=lambda x: (x[1]["primer_reverse_confidence"], x[1]["mean_error_rate"]), reverse=True)
        self.ranking_reverse = [r[0] for r in ranking_reverse]
        self.run_infos["ranking_computed"] == True
    
    def get_best_sequences(self, n=100, key="forward", skipped=0, truncate_sequences=True):
        """
        Get the `n` best sequences according to the [forward/reverse] primer
        :param n: number of sequences [int]
        :param key: "forward" or "reverse" [str]
        :param skipped: number of best sequences "skipped" before getting n sequences [int]
        :return: n best sequences according to the [forward/reverse] primer
        """
        if self.run_infos["ranking_computed"] == False:
            self.compute_sequence_ranking()
            warnings.warn("[MultiSequence.get_best_sequences] compute_sequence_ranking() has been called to get the best sequences")
        if key == "forward":
            ranking = self.ranking_forward
        elif key == "reverse":
            ranking = self.ranking_reverse
        else:
            warnings.warn(f"[MultiSequence.get_best_sequences] 'key' argument must be equal to either 'forward' or 'reverse'")
        indices = ranking[skipped:skipped+n]
        best_sequences = copy.deepcopy([self.sequences[index] for index in indices])
        if truncate_sequences:
            for Seq_id in range(len(best_sequences)):
                best_sequences[Seq_id].truncate_sequence_using_best_primer_alignment()
        BestSequences = MultiSequence()
        BestSequences.sequences = best_sequences
        return BestSequences
    
    def apply_threshold_error_rate_max(self, max_error_rate=0.06):
        """
        Apply a threshold (maximum allowed for the mean_error_rate) and remove the Sequence() not in the allowed range
        :param max_error_rate: maximum mean_error_rate allowed
        :return: None
        """
        self.sequences = [seq for seq in self.sequences if seq.get_scores()["mean_error_rate"] <= max_error_rate]
        self.compute_sequence_ranking()

    def apply_threshold_primer_alignment_confidence_min(self, min_primer_alignment_confidence=0.8):
        self.sequences = [seq for seq in self.sequences if seq.get_scores()["primer_alignment_confidence"] >= min_primer_alignment_confidence]
        self.compute_sequence_ranking()
    
    def print_stats_sequences(self):
        """
        Print informations/stats about the sequences (base composition, sequence lengths)
        :return: None
        """
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        seq_lengths = []
        base_A = 0
        base_T = 0
        base_G = 0
        base_C = 0
        for seq in self.sequences:
            seq_lengths.append(len(seq.sequence))
            base_A += seq.sequence.count('A')
            base_T += seq.sequence.count('T')
            base_G += seq.sequence.count('G')
            base_C += seq.sequence.count('C')
        base_total = base_A + base_T + base_G + base_C
        seq_mean_length = np.mean(np.array(seq_lengths))
        seq_std_length = np.std(np.array(seq_lengths))
        seq_median_length = np.median(np.array(seq_lengths))
        print(f"Stats about the sequences (n={len(self.sequences)}):")
        print(f"\tPercentage of A: {round(100 * base_A / base_total, 2)}%")
        print(f"\tPercentage of T: {round(100 * base_T / base_total, 2)}%")
        print(f"\tPercentage of G: {round(100 * base_G / base_total, 2)}%")
        print(f"\tPercentage of C: {round(100 * base_C / base_total, 2)}%")
        print(f"\tMean length:   {round(seq_mean_length, 1)}")
        print(f"\tStd length:    {round(seq_std_length, 1)}")
        print(f"\tMedian length: {round(seq_median_length, 1)}")
        print(f"\tPercentiles length: ", end='')
        for percent in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
            print(f"{round(np.percentile(seq_lengths, percent))} ", end='')
        print(f"\n\t{20 * ' '}1%  5%  10% 25% 50% 75% 90% 95% 99%")
    
    def print_stats_error_rates(self):
        """
        Print informations/stats about the sequence error rates (sequencing error rates)
        :return: None
        """
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        seq_error_rates = []
        for seq in self.sequences:
            mean_error_rate = seq.get_scores()["mean_error_rate"]
            if mean_error_rate is not None: seq_error_rates.append(mean_error_rate)
        if len(seq_error_rates) == 0:
            print("No information about the phred scores")
            return None
        count_error_rate = np.round(100 * np.array(seq_error_rates))
        values, counts = np.unique(count_error_rate, return_counts=True)
        cumulative_percentage = 0.0
        print(f"Stats about the sequencing errors (n={len(self.sequences)}):")
        print(f"\tInfos: number of sequences (proportion|cumulative)")
        for i in range(len(values)):
            percentage = 100 * counts[i] / count_error_rate.size
            cumulative_percentage += percentage
            print(f"\t{round(values[i])}% of errors: {counts[i]} ({int(percentage)}%|{int(cumulative_percentage)}%)")
        error_rate_median = np.median(np.array(seq_error_rates))
        error_rate_mean = np.mean(np.array(seq_error_rates))
        error_rate_std = np.std(np.array(seq_error_rates))
        print(f"\tMedian error rate: {round(100 * error_rate_median, 2)}%")
        print(f"\tMean error rate: {round(100 * error_rate_mean, 2)}%")
        print(f"\tStd error rate: {round(100 * error_rate_std, 2)}%")
    
    def print_stats(self):
        """
        Print all implemented informations/stats about MultiSequence() instance
        :return: None
        """
        if len(self.sequences) == 0:
            warnings.warn("No sequence found in MultiSequence() instance.")
            return None
        self.print_stats_sequences()
        self.print_stats_error_rates()
    
    def __str__(self):
        """
        print(MultiSequence()) uses this function to format the MultiSequence() into string representation.
        When printing a MultiSequence() instance: print_stats() is called
        :return: empty string [str]
        """
        print("/!\ Printing a MultiSequence() instance calls the method print_stats().")
        self.print_stats()
        return ""