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

    def tearDown(self):
        """ Remove all sequences from the dictionary """

    ###########################################################################
    #   T E S T  C A S E S
    ###########################################################################

    def test_file_open(self):
        """ Test file open, read, dictionary"""
        name = "gi|142022655|gb|EQ086233.1|43"
        self.assertEqual(name in self.test_seqs, True)
        fs = sequences.FastaSeq()
        fs.buildDict()
        self.assertEqual(fs.getSeq(name), self.test_seqs[name])

    def test_number_recs(self):
        """ It should return the number of records """
        num_recs = len(self.test_seqs)
        fs = sequences.FastaSeq()
        fs.buildDict()
        self.assertEqual(num_recs, fs.numRecords())

    def test_seq_length(self):
        """ It should return the length of a record's sequence """
        name = "gi|142022655|gb|EQ086233.1|43"
        fs = sequences.FastaSeq()
        fs.buildDict()
        self.assertEqual(len(self.test_seqs[name]),fs.getLength(name))

    def test_all_lengths(self):
        """ It should return the length of all sequences """
        fs = sequences.FastaSeq()
        fs.buildDict()
        lengths = [len(x) for x in list(self.test_seqs.values())]
        self.assertEqual(sorted(lengths), sorted(fs.getAllLengths()))

    def test_longest_seqs(self):
        """ It should return a dict Length:[identifiers] """
        fs = sequences.FastaSeq()
        fs.buildDict()
        max_length = 4805
        longest = {max_length:[]}
        for l in iter(self.test_seqs):
            if len(self.test_seqs[l]) == max_length:
                longest[max_length] += [l]
        self.assertEqual(longest[max_length], fs.getLongest()[max_length])
        


        