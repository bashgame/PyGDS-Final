"""
Test Cases for DNA sequences
"""
from unittest import TestCase
from source import sequences


class TestSequences(TestCase):
    """ Tests for sequences.py """
    FILENAME = 'fixtures/dna.example.fasta'
    seqs = {}

    @classmethod
    def setUpClass(cls):
        """ Open cls.FILENAME """
        cls.file = open(cls.FILENAME,r)

    @classmethod
    def tearDownClass(cls):
        """ Close file """
        cls.file.close()

    def setUp(self):
        """ Set up our dictionary of sequences """
        for line in self.file.readline():
            if line.startswith('>'):
                name = line[1:]
                self.seqs[name] = ''
            else:
                self.seqs[name] = self.seqs[name] + line.strip()

    def tearDown(self):
        """ Remove all sequences from the dictionary """
        self.seqs.clear()

    ###########################################################################
    #   T E S T  C A S E S
    ###########################################################################

    def test_file_open(self):
        """ Test file open, read, dictionary"""
        first_name = self.seqs.list()[0]
        self.assertEqual(first_name, "gi|142022655|gb|EQ086233.1|43 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence")
