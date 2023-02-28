# from Bio import SeqIO
# from Bio import Seq

# seq_recs = {record.id: record for record in SeqIO.parse(FILENAME, 'fasta')}


class FastaSeq():
    sequences = {}

    @classmethod
    def buildDict(self, filename):
        """ buildDict builds a dictionary of name: sequence pairs given a fasta formatted file

        Args:
            filename (str): the name of the fasta file to open
        """
        self.sequences.clear()
        try:
            file = open(filename)
        except FileNotFoundError:
            print(f"File {filename} not found!")
            return
        else:
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
        """ numRecords() is a function to get the number of records in the class dictionary

        Returns:
            int num_records: the number of records in the class dictionary returned by len(dict)
        """
        return len(self.sequences)

    @classmethod
    def getLength(self, name):
        """ getLength returns the length of a sequence associated with a name

        Args:
            name (str): The name of the dictionary key

        Returns:
            int length: the length of the sequence associated with name in the class dictionary
        """
        return len(self.sequences[name])

    @classmethod
    def getAllLengths(self):
        """ getAllLengths returns a list of the lengths of all sequences in the class dictionary

        Returns:
            [ int lengths ]: a list of ints representing the lengths of all sequences
        """
        return [len(x) for x in list(self.sequences.values())]

    @classmethod
    def getLongest(self):
        """ getLongest returns a dictionary with one key, the max length of sequences in the dictionary,
        and a value consisting of the list of names that have a sequence of that length

        Returns:
            { int max_length: [ str names ] }: The maximum length of all sequences in the dictionary:
            A list of all names with a sequence of that length
        """
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
        """ getShortest gets the length of the shortest sequence and a list of all names with a sequence of that length

        Returns:
            { int min_length: [ str names ] }: the length of the shortest sequence in the dictionary:
            a list of all names with sequences of that length
        """
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
    def getStopCodons(self, sequence):
        """ getStopCodons looks for the subsequences 'tga', 'tag', and 'taa' in sequence. When it finds one of these
        it locates which reading frame the codon is located in and stores the index in a list associated with that
        reading frame through a dictionary. The reading frames are numbered 0-2 rather than 1-3.

        Args:
            sequence (str): The sequence of nucleotides to search for stop codons

        Returns:
            { int reading_frame: [ int index ]}: A dictionary with keys consisting of reading frame 0, 1, 2 and values
            consisting of lists of indices of stop codons on those reading frames.
        """
        stop_codons = ['tga', 'tag', 'taa']
        stops_dict = {x: [] for x in range(3)}
        seq = sequence.lower()
        for idx in range(len(seq)):             # iterate over each codon in the sequence
            codon = seq[idx:idx+3]
            if codon in stop_codons:
                stops_dict[idx % 3] += [idx]    # idx % 3 gives the frame, add idx to the frame's list

        return stops_dict

    @classmethod
    def getStartCodons(self, sequence):
        """ getStartCodons looks for the subsequence 'atg' in sequence. When it finds this it locates which reading
        frame the codon is located in and stores the index in a list associated with that reading frame through a
        dictionary. The reading frames are numbered 0-2 rather than 1-3.

        Args:
            sequence (str): A sequence of nucleotides to search for start codons

        Returns:
            { int reading_frame: [ int index ]}: A dictionary with keys consisting of reading frame 0, 1, 2 and values
            consisting of lists of indices of start codons on those reading frames.
        """
        starts_dict = {x: [] for x in range(3)}
        seq = sequence.lower()
        for idx in range(len(seq)):
            codon = seq[idx:idx+3]
            if codon == 'atg':
                starts_dict[idx % 3] += [idx]

        return starts_dict

    @classmethod
    def getLongestORF(self, sequence):
        """ getLongestORF calls getStartCodons and getStopCodons and then computes the longest possible Open Reading
        Frame by computing the difference between the first start codon and the last stop codon on each reading frame.
        An Open Reading Frame must begin with a start codon and end with a stop codon on a reading frame.

        Args:
            sequence (str): A string of nucleotides to search for Open Reading Frames.

        Returns:
            { 'length': int longest, 'index': int index}: A dictionary containing the length and index of the longest
            ORF in the sequence
        """
        starts = self.getStartCodons(sequence)
        stops = self.getStopCodons(sequence)
        lngst = 0
        for frame in iter(starts):
            if stops[frame] == [] or starts[frame] == []:
                continue
            length = stops[frame][-1] - starts[frame][0]
            if length > lngst:
                lngst = length
                lng_idx = starts[frame][0]
        orf_dict = {'length': lngst, 'index': lng_idx}
        return orf_dict

    @classmethod
    def getFileLongestORF(self):
        """ getFileLongestORF calls getLongestORF for each sequence in the dictionary. It then stores the longest ORF
        in each sequence in a new dictionary, orf_dict, with the key being the name from the class dictionary and the
        corresponding value being a length and index for the longest ORF in that name's sequence. The method also tracks
        the length and index of the longest ORF of all the sequences, as well as the name associated with that sequence.
        It then returns all of this to the caller using a dictionary with 'name', 'length', 'index', and 'data' keys.
        The 'name', 'length' and 'index' keys all correspond to the longest ORF in the class dictionary. The 'data' key
        corresponds to the dictionary of longest ORFs for each name and sequence in the class dictionary.

        Returns:
            { 'name': str lgst_name, 'length': int longest, 'index': lgst_idx, 'data': { str name { 'length': int length,
            'index': int index} } }: A dictionary consisting of the name, length, and starting index of the longest ORF
            in the class dictionary, as well as a separate dictionary containing the name, length and starting incex of
            the longest ORF for each sequence in the class dictionary object
        """
        longest = 0
        lgst_name = ''
        orf_dict = {}
        for name, seq in self.sequences.items():
            result = self.getLongestORF(seq)
            orf_dict[name] = result
            if result["length"] > longest:
                longest = result["length"]
                lgst_name = name
                lgst_idx = result["index"]
        ret_dict = {"name": lgst_name, "length": longest, "index": lgst_idx, "data": orf_dict}
        return ret_dict

    @classmethod
    def getRepeats(self, sequence, length):
        """ getRepeats searches for repeat sequences of length in sequence

        Args:
            sequence (str): A sequence of nucleotides to search
            length (int): The length of subsequences to search for

        Returns:
            { str substr: int repeats }: A dictionary consisting of keys substr and values repeats. The dictionary
            will contain a key for each substring of length in sequence, and the value will contain the number of
            times that substring has been repeated. Note that if a substring only appears once in the sequence it
            will have a value of 0 in the dictionary.
        """
        repeats = {}
        for idx in range(len(sequence)):
            substr = sequence[idx:idx+length].lower()
            if substr in list(repeats):
                repeats[substr] += 1
            elif len(substr) == length:
                repeats[substr] = 0                 # We are counting repeats of the substring
        return repeats

    @classmethod
    def getMostRepeats(self, rep_dict):
        """ getMostRepeats searches rep_dict for the most common repeat in the dictionary

        Args:
            rep_dict ( { str: int } ): A dictionary with keys substring and int repetitions

        Returns:
            { str most_common: int most_reps }: This is a single key, value pair with the
            most common string in rep_dict as key and the value being the same value
            associated with that key in rep_dict.
        """
        most_common = ''
        most_reps = 0
        for key, val in rep_dict.items():
            if val > most_reps:
                most_reps = val
                most_common = key
        return {most_common: most_reps}

    @classmethod
    def getMultiSeqRepeats(self, seq_dict, length):
        """ getMultiSeqRepeats repeatedly calls getRepeats for each sequence in seq_dict, searching for repeat
        substrings of length. It then combines the dictionaries that getRepeats returns into one dictionary
        containing the totals for all substrings of length found in each sequence in seq_dict

        Args:
            seq_dict ({ str name: str sequence}): A dictionary of name: sequence pairs
            length (int): The length of subsequences to search for

        Returns:
            { str substr: int repeats }: A dictionary consisting of keys substr and values repeats. The dictionary
            will contain a key for each substring of length in all sequences in seq_dict, and the value will contain
            the number of times that substring has been repeated. Note that if a substring only appears once in all
            sequences it will have a value of 0 in the dictionary.
        """
        reps_dict = {}
        for seq in seq_dict.values():
            tmp_dict = self.getRepeats(seq, length)
            for key, val in tmp_dict.items():
                if key in reps_dict:
                    reps_dict[key] += val
                else:
                    reps_dict[key] = val
        return reps_dict


def main():
    return
