import Bio

FILENAME = 'tests/fixtures/dna.example.fasta'

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
        longest = {max_length:[]}
        counts = 0
        for n in iter(self.sequences):
            if len(self.sequences[n]) == max_length:
                longest[max_length] += [n]
                counts += 1
                if counts == num_max:
                    break
        return longest

    @classmethod
    def getSeq(self, name):
        return self.sequences[name]

def main():
    return

