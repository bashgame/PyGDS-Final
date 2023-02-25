# from Bio import SeqIO
#from Bio import Seq

FILENAME = 'tests/fixtures/dna.example.fasta'
# seq_recs = {record.id: record for record in SeqIO.parse(FILENAME, 'fasta')}


class FastaSeq():
    sequences = {}

    @classmethod
    def buildDict(self):
        self.sequences.clear()
        file = open(FILENAME)
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):
                words = line.split()
                name = words[0][1:]
                self.sequences[name] = ''
            else:
                self.sequences[name] = self.sequences[name] + line
        file.close()

    @classmethod
    def numRecords(self):
        return len(self.sequences)

    @classmethod
    def getLength(self, name):
        return len(self.sequences[name])

    @classmethod
    def getAllLengths(self):
        return [len(x) for x in list(self.sequences.values())]

    @classmethod
    def getLongest(self):
        lengths = self.getAllLengths()
        max_length = max(lengths)
        num_max = lengths.count(max_length)
        longest = {max_length: []}
        counts = 0
        for n in iter(self.sequences):
            if len(self.sequences[n]) == max_length:
                longest[max_length] += [n]
                counts += 1
                if counts == num_max:
                    break
        return longest

    @classmethod
    def getShortest(self):
        lengths = self.getAllLengths()
        min_length = min(lengths)
        shortest = []
        for n in iter(self.sequences):
            if len(self.sequences[n]) == min_length:
                shortest += [n]

        shortest = {min_length: shortest}
        return shortest

    @classmethod
    def getSeq(self, name):
        return self.sequences[name]

    @classmethod
    def getStopCodons(self, name):
        stop_codons = ['tga','tag','taa']
        stops_dict = {x: [] for x in range(3)}
        seq = self.sequences[name].lower()
        for idx in range(len(seq)): # iterate over each codon in the sequence
            codon = seq[idx:idx+3]
            if codon in stop_codons:
                stops_dict[idx % 3] += [idx]  # idx % 3 gives the frame, add idx to the frame's list

        return stops_dict

    @classmethod
    def getStartCodons(self, name):
        starts_dict = {x: [] for x in range(3)}
        seq = self.sequences[name].lower()
        for idx in range(len(seq)):
            codon = seq[idx:idx+3]
            if codon == 'atg':
                starts_dict[idx % 3] += [idx]

        return starts_dict

def main():
    return
