"""
Test Cases for DNA sequences
"""
from unittest import TestCase
from source import sequences


class TestSequences(TestCase):
    """ Tests for sequences.py """
    FILENAME = 'tests/fixtures/dna.example.fasta'
    seqs = {}
    test_seqs = {}
    fs = sequences.FastaSeq()

    @classmethod
    def setUpClass(cls):
        """ Open cls.FILENAME """
        cls.file = open(cls.FILENAME,'r')
        for line in cls.file:
            line = line.rstrip()
            if line.startswith('>'):
                words = line.split()
                name = words[0][1:]
                cls.seqs[name] = ''
            else:
                cls.seqs[name] = cls.seqs[name] + line
        cls.file.close()
    
    @classmethod
    def tearDownClass(cls):
        """ Close file """
        cls.seqs.clear()

    def setUp(self):
        """ Set up our dictionary of sequences """
        self.test_seqs.clear()
        self.test_seqs = self.seqs
        self.fs.buildDict(self.FILENAME)

    def tearDown(self):
        """ Remove all sequences from the dictionary """

    ###########################################################################
    #   T E S T  C A S E S
    ###########################################################################

    def test_file_open(self):
        """ Test file open, read, dictionary"""
        name = "gi|142022655|gb|EQ086233.1|43"
        self.assertEqual(name in self.test_seqs, True)
        self.assertEqual(self.fs.getSeq(name), self.test_seqs[name])

    def test_number_recs(self):
        """ It should return the number of records """
        num_recs = len(self.test_seqs)
        self.assertEqual(num_recs, self.fs.numRecords())

    def test_seq_length(self):
        """ It should return the length of a record's sequence """
        name = "gi|142022655|gb|EQ086233.1|43"
        self.assertEqual(len(self.test_seqs[name]), self.fs.getLength(name))

    def test_all_lengths(self):
        """ It should return the length of all sequences """
        lengths = [len(x) for x in list(self.test_seqs.values())]
        self.assertEqual(sorted(lengths), sorted(self.fs.getAllLengths()))

    def test_longest_seqs(self):
        """ It should return a dict max_length:[identifiers] """
        max_length = 4805
        longest = {max_length:[]}
        for l in iter(self.test_seqs):
            if len(self.test_seqs[l]) == max_length:
                longest[max_length] += [l]
        self.assertEqual(longest[max_length], self.fs.getLongest()[max_length])
        
    def test_shortest_seqs(self):
        """ It should return a dict min_length:[identifiers] """
        min_length = 512
        shortest = []
        for l in iter(self.test_seqs):
            if len(self.test_seqs[l]) == min_length:
                shortest += [l]
        self.assertEqual(shortest, self.fs.getShortest()[min_length])

    def test_stop_codons(self):
        """ Given an id, it should return a dict with a list of indices for stop codons for each frame """
        name = "gi|142022655|gb|EQ086233.1|43"
        stop_codons = ['tga','tag','taa']
        stops = 0
        for stc in stop_codons:
            stops += self.test_seqs[name].lower().count(stc)
        stc_seq = self.fs.getSeq(name)
        stc_ret = self.fs.getStopCodons(stc_seq)
        self.assertEqual(stops, len(stc_ret[0]) + len(stc_ret[1]) + len(stc_ret[2]))

    def test_start_codons(self):
        """ Given an id, it should return a dict with a list of indices for start codons for each frame """
        name = "gi|142022655|gb|EQ086233.1|43"
        starts = self.test_seqs[name].lower().count('atg')
        stc_seq = self.fs.getSeq(name)
        stc_ret = self.fs.getStartCodons(stc_seq)
        self.assertEqual(starts, len(stc_ret[0]) + len(stc_ret[1]) + len(stc_ret[2]))

    def test_longest_orf_in_seq(self):
        """ Given an id, return the longest ORF on that sequence """
        name = "gi|142022655|gb|EQ086233.1|43"
        sequence = self.fs.getSeq(name)
        starts = self.fs.getStartCodons(sequence)
        stops = self.fs.getStopCodons(sequence)
        orf_dict = {}
        for frame in range(3):
            lngst = 0
            lng_idx = -1
            idx1 = 0
            if stops[frame] != [] and starts[frame] != []:
                for idx2 in range(len(stops[frame])):
                    if idx1 < len(starts[frame]) and starts[frame][idx1] < stops[frame][idx2]:
                        length = 3 + stops[frame][idx2] - starts[frame][idx1]
                        if length > lngst:
                            lngst = length
                            lng_idx = starts[frame][idx1]
                        while idx1 < len(starts[frame]) and starts[frame][idx1] < stops[frame][idx2]:
                            idx1 += 1
            orf_dict[frame] = {'length': lngst, 'index': lng_idx}
        max_len = orf_dict[0]['length']
        max_idx = orf_dict[0]['index']
        result = self.fs.getLongestORF(sequence)
        self.assertEqual(max_len, result[0]["length"])
        self.assertEqual(max_idx, result[0]["index"])

    def test_longest_orf_in_file(self):
        """ Return id, length, and index of the longest ORF in the file """
        orf_dict = {}
        ret_dict = {}
        longest = [0, 0, 0]
        lgst_name = ['', '', '']
        lgst_idx = [0, 0, 0]
        for name, seq in self.fs.sequences.items():
            result = self.fs.getLongestORF(seq)
            orf_dict[name] = result
            for frame in range(3):
                if result[frame]["length"] > longest[frame]:
                    longest[frame] = result[frame]["length"]
                    lgst_name[frame] = name
                    lgst_idx[frame] = result[frame]["index"]
                ret_dict[frame] = {"name": lgst_name[frame], "length": longest[frame], "position": lgst_idx[frame] + 1}
        max_name = ret_dict[0]['name']
        max_len = ret_dict[0]['length']
        max_idx = ret_dict[0]['position']
        result = self.fs.getFileLongestORF()
        self.assertEqual(max_name, result[0]["name"])
        self.assertEqual(max_len, result[0]["length"])
        self.assertEqual(max_idx, result[0]["position"])

    def test_repeats(self):
        """ Given a sequence and length n, it should return all repeats
        of length n and the number of their occurrences """ 
        name = "gi|142022655|gb|EQ086233.1|43"
        test_str = 'ACACAGGGACACA'
        ac_reps = 4
        ca_reps = 4
        gg_reps = 2
        aca_reps = 4
        cac_reps = 2
        acac_reps = 2
        caca_reps = 2
        acaca_reps = 2
        result = self.fs.getRepeats(test_str, 2)
        self.assertEqual(ac_reps, result["ac"])
        result = self.fs.getRepeats(test_str, 3)
        self.assertEqual(aca_reps, result["aca"])
        self.assertEqual(cac_reps, result["cac"])
        result = self.fs.getRepeats(self.test_seqs[name], 3)
        self.assertNotEqual(0, result["aca"])

    def test_get_frequent_repeat(self):
        """ Given a dictionary of repeats and their occurrences, it should return
        the repeat with the highest occurrence  """
        test_str = 'ACACAGGGACACA'
        test_reps = self.fs.getRepeats(test_str, 2)
        exp_result = {"ac": 4}
        result = self.fs.getMostRepeats(test_reps)
        self.assertEqual(exp_result, result)

    def test_get_multiseq_repeats(self):
        """ Given a dictionary of sequences and a length n, it should return the
        repeats of substrings of that length across all sequences """
        test_str = 'ACACAGGGACACA'
        aca_reps = 4
        exp_result = 4 * aca_reps
        test_names = ["Alice", "Bob", "Cho", "Darwish"]
        test_dict = {x: test_str for x in test_names}
        result = self.fs.getMultiSeqRepeats(test_dict, 3)
        self.assertEqual(exp_result, result['aca'])

    def test_get_file_repeats(self):
        """ It should return the most common repeat in the file """
        result = self.fs.getMultiSeqRepeats(self.fs.sequences, 3)
        result = self.fs.getMostRepeats(result)
        self.assertNotEqual(0, result)
            