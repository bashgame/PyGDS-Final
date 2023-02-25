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
        stc_ret = self.fs.getStopCodons(name)
        self.assertEqual(stops, len(stc_ret[0]) + len(stc_ret[1]) + len(stc_ret[2]))

    def test_start_codons(self):
        """ Given an id, it should return a dict with a list of indices for start codons for each frame """
        name = "gi|142022655|gb|EQ086233.1|43"
        starts = self.test_seqs[name].lower().count('atg')
        stc_ret = self.fs.getStartCodons(name)
        self.assertEqual(starts, len(stc_ret[0]) + len(stc_ret[1]) + len(stc_ret[2]))

    def test_longest_orf_in_seq(self):
        """ Given an id, return the longest ORF on that sequence """
        name = "gi|142022655|gb|EQ086233.1|43"
        start_codon = 'atg'
        stop_codons = ['tga','tag','taa']
        start_fx = [0, 0, 0]
        stop_fx = [0, 0, 0]
        seq = self.seqs[name].lower()
        for frame in range(3):
            for idx in range(frame,len(seq),3):
                codon = seq[idx:idx+3]
                if codon == start_codon:
                    start_fx[frame] = idx
                    break        
        end = len(seq)
        while 0 in stop_fx:
            idx = max(seq.rfind(stop_codons[0], 0, end), seq.rfind(stop_codons[1], 0, end), seq.rfind(stop_codons[2], 0, end))
            if idx > stop_fx[idx % 3]:
                stop_fx[idx % 3] = idx
            end = idx - 1
        max_len = 0
        max_idx = 0
        for frame in range(3):
            length = stop_fx[frame] - start_fx[frame]
            if length > max_len:
                max_len = length
                max_idx = start_fx[frame]
        result = self.fs.getLongestORF(name)
        self.assertEqual(max_len, result["length"])
        self.assertEqual(max_idx, result["index"])

    def test_longest_orf_in_file(self):
        """ Return id, length, and index of the longest ORF in the file """
        start_codon = 'atg'
        stop_codons = ['tga','tag','taa']
        glob_max_len = 0
        glob_max_idx = 0
        glob_max_name = ''
        for name in iter(self.test_seqs):
            start_fx = [0, 0, 0]
            stop_fx = [0, 0, 0]
            seq = self.seqs[name].lower()
            for frame in range(3):
                for idx in range(frame,len(seq),3):
                    codon = seq[idx:idx+3]
                    if codon == start_codon:
                        start_fx[frame] = idx
                        break        
            end = len(seq)
            while 0 in stop_fx and end > 0:
                idx = max(seq.rfind(stop_codons[0], 0, end), seq.rfind(stop_codons[1], 0, end), seq.rfind(stop_codons[2], 0, end))
                if idx > stop_fx[idx % 3]:
                    stop_fx[idx % 3] = idx
                end = idx - 1
            loc_max_len = 0
            loc_max_idx = 0
            for frame in range(3):
                length = stop_fx[frame] - start_fx[frame]
                if length > loc_max_len:
                    loc_max_len = length
                    loc_max_idx = start_fx[frame]
            if loc_max_len > glob_max_len:
                glob_max_len = loc_max_len
                glob_max_name = name
                glob_max_idx = loc_max_idx
        result = self.fs.getFileLongestORF()
        self.assertEqual(glob_max_name, result["name"])
        self.assertEqual(glob_max_len, result["length"])
        self.assertEqual(glob_max_idx, result["index"])

