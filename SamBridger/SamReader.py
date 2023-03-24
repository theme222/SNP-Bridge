import pysam
from Bio import SeqIO
from random import randint, choice
from Bio.Align import PairwiseAligner
from TechnicalTools import multi_delete, log
import numpy as np
from collections import Counter


class SNP:
    """ class to store SNP data type containing allele and global position"""

    def __init__(self, global_start_pos, local_start_pos=None):
        self.allele = []
        self.global_start_pos = global_start_pos
        if local_start_pos is None:
            self.local_start_pos = global_start_pos
        else:
            self.local_start_pos = local_start_pos

    def __repr__(self):
        return f'{self.global_start_pos}'

    '''
    def get_allele(self, *args):
        if not self.allele:  # self.allele == []
            for i in args:
                self.allele.append(i.read[self.start_pos])
        return self.allele
    '''

    def snp_checker(self, snplist):  # check whether a snp is in a lsist
        try:
            index = [i.global_start_pos for i in snplist].index(self.global_start_pos)  # index within snplist
            return True, index
        except ValueError:
            return False, None

    @staticmethod
    def key(snp):
        return snp.global_start_pos


class DNA:
    """
    class to store DNA type information
    """

    nucleotides = ['A', 'T', 'G', 'C']
    current_insertid = 1  # keeping track and making sure insert Ids aren't the same
    global_snp = []  # all snps in the current DNA strand

    def __init__(self, read=np.array([], dtype=str), snp=None, start_pos=0, insertions=None, insertion_pos=None):
        if snp is None: snp = []
        if insertion_pos is None:  insertion_pos = np.array([], dtype=int)
        if insertions is None:  insertions = np.array([], dtype=str)
        self.read = read  # Numpy array
        self.snp = snp  # Normal Python list containing SNP objects?
        self.start_pos = start_pos  # Int
        self.insertid = 0  # Int
        self.insertions = insertions  # list of Char
        self.insertion_pos = insertion_pos  # list of numbers

    def __repr__(self):
        return f'"{self.start_pos}"'
        #  return f'{list(self.read)} \n mutates at {self.snp}'

    def simulate(self, size):
        self.read = np.random.choice(DNA.nucleotides, size + 1)

    @staticmethod
    def mutate(type_a, type_b, snp_count):
        snp_pos = []
        for i in range(snp_count):
            r = randint(1, type_a.size) - 1
            if r not in snp_pos:
                snp_pos.append(r)

        for l in snp_pos:
            current_type = choice([type_a, type_b])
            # making sure the choice isn't the same as the original
            DNA.nucleotides.remove(current_type.read[l])
            current_type.read[l] = choice(DNA.nucleotides)
            DNA.nucleotides = ['A', 'T', 'G', 'C']

            current_type.snp.append(SNP(l))
            DNA.global_snp.append(SNP(l))
        type_a.snp.sort(key=SNP.key)
        type_b.snp.sort(key=SNP.key)
        DNA.global_snp.sort(key=SNP.key)

    """
    Adenine == 1
    Thymine == 2
    Cytosine == 3
    Guanine == 4
    """

    def to_int(self):

        if self.read.dtype == int:
            return self.read
        elif self.read.dtype == '<U1':

            if len(self.read.shape) == 1:
                _temp = np.array([], dtype=int)  # return variable
                for acid in self.read:
                    if acid == 'A' or acid == 'a':
                        _temp = np.append(_temp, 1)
                    elif acid == 'T' or acid == 't':
                        _temp = np.append(_temp, 2)
                    elif acid == 'C' or acid == 'c':
                        _temp = np.append(_temp, 3)
                    elif acid == 'G' or acid == 'g':
                        _temp = np.append(_temp, 4)
                    else:
                        log('error', 'Unknown character in to int function', acid)
                return _temp
            else:
                _temp2 = np.zeros((0, self.read.shape[1]), dtype=int)
                for read in self.read:
                    _temp = np.array([], dtype=int)  # return variable
                    for acid in read:
                        if acid == 'A' or acid == 'a':
                            _temp = np.append(_temp, 1)
                        elif acid == 'T' or acid == 't':
                            _temp = np.append(_temp, 2)
                        elif acid == 'C' or acid == 'c':
                            _temp = np.append(_temp, 3)
                        elif acid == 'G' or acid == 'g':
                            _temp = np.append(_temp, 4)
                        else:
                            log('error', 'Unknown character in to int function', acid)
                    _temp2 = np.concatenate((_temp2, [_temp]))
                return _temp2

    def to_str(self):

        if self.read.dtype == '<U1':
            return self.read
        elif self.read.dtype == int:

            if len(self.read.shape) == 1:
                _temp = np.array([], dtype=int)  # return variable
                for acid in self.read:
                    if acid == 1:
                        _temp = np.append(_temp, 'A')
                    elif acid == 2:
                        _temp = np.append(_temp, 'T')
                    elif acid == 3:
                        _temp = np.append(_temp, 'C')
                    elif acid == 4:
                        _temp = np.append(_temp, 'G')
                    else:
                        log('error', 'Unknown character in to asc function', acid)
                return _temp
            else:
                _temp2 = np.zeros((0, self.read.shape[1]), dtype=int)
                for read in self.read:
                    _temp = np.array([], dtype=int)  # return variable
                    for acid in read:
                        if acid == 1:
                            _temp = np.append(_temp, 'A')
                        elif acid == 2:
                            _temp = np.append(_temp, 'T')
                        elif acid == 3:
                            _temp = np.append(_temp, 'C')
                        elif acid == 4:
                            _temp = np.append(_temp, 'G')
                        else:
                            log('error', 'Unknown character in to asc function', acid)
                    _temp2 = np.concatenate((_temp2, [_temp]))
                return _temp2

    @property
    def size(self):
        return len(self.read)

    @classmethod
    def derive(cls, obj):
        return cls(read=obj.read.copy())

    @classmethod
    def obtain_read(cls, obj, read_size, insert_size, index=None, insertid=0):  # Makes one insert
        if index is None:
            index = randint(0, obj.size - insert_size)

        # Goes through the global snp list and finds a SNP object to then make a new object with global and local position
        forwardSNP = [SNP(i.global_start_pos, local_start_pos=i.global_start_pos - index) for i in obj.global_snp
                      if index < i.global_start_pos < index + read_size]

        reverseSNP = [SNP(i.global_start_pos, local_start_pos=i.global_start_pos - index - insert_size + read_size) for
                      i in
                      obj.global_snp
                      if index + insert_size - read_size < i.global_start_pos < index + insert_size]

        forward = cls(read=obj.read[index:index + read_size], snp=forwardSNP, start_pos=index)
        reverse = cls(read=obj.read[index + insert_size - read_size:index + insert_size], snp=reverseSNP,
                      start_pos=index + insert_size - read_size)
        forward.insertid = insertid
        reverse.insertid = insertid
        return [forward, reverse]

    @classmethod
    def snp_purger(cls, read_list):  # finding snps that aren't snps

        #  Purge useless snps
        snp_values = []
        final_snp_list = []
        for index, snp in enumerate(cls.global_snp):
            for read in read_list:
                if read.start_pos <= snp.global_start_pos <= read.start_pos + read.size - 1:
                    snp_values.append(read.read[snp.global_start_pos - read.start_pos])

            # Combine all words together
            chars = "".join(snp_values)
            # Count all characters
            counter = Counter(chars)
            # Get frequency and put in list
            total_occ = [(char, occ) for char, occ in counter.most_common()]
            if len(total_occ) == 2:  # only include a snp when it has two values
                if total_occ[1][1] / total_occ[0][1] > 0.4:
                    final_snp_list.append(snp)
            snp_values = []

        cls.global_snp = final_snp_list
        cls.global_snp.sort(key=SNP.key)
        log('critical', cls.global_snp)

    @staticmethod
    def blend(strand, avg_depth, read_size, insert_size):  # Makes multiple inserts
        return_list = np.array([], dtype=object)
        insert_start_pos = np.random.randint(0, strand.size - insert_size,
                                             round((strand.size * avg_depth) / read_size / 2))
        insert_start_pos.sort()
        for i in insert_start_pos:
            return_list = np.append(return_list, DNA.obtain_read(
                strand, read_size, insert_size, index=i, insertid=DNA.current_insertid))
            DNA.current_insertid += 1
        return return_list

    @staticmethod
    def find_rev_by_insertid(read, shuffled_reads):
        insertid = read.insertid
        for i in shuffled_reads:
            if i.insertid == insertid and i != read:
                return i


def sam_aligner(reference, sam_list, sam_start_pos):
    """
    reference: String of reference gene
    sam_list: list of strings
    sam_start_pos: list of integers of starting position 0
    return : DNA object with SNPs inside
    """

    # return variables
    snp_compare_list = []
    snp_list = []
    mod_sam_list = np.array([], dtype=object)

    genewithposlist = []  # a list of lists with [sequence, starting position]
    if len(sam_list) == len(sam_start_pos):
        for i in range(len(sam_list)):
            genewithposlist.append([sam_list[i], sam_start_pos[i]])

    for sam in genewithposlist:
        start_pos = sam[1]
        sam_seq = sam[0]

        # making reference that is slightly bigger that target sequence so it can account for ins/del
        seq1 = reference[start_pos:start_pos + len(sam_seq) + 20]

        seq2 = sam_seq  # sam sequence

        # =========== AI code incoming ===========
        # Create a PairwiseAligner instance
        aligner = PairwiseAligner()

        # Set the mode to global
        aligner.mode = 'global'

        # Set the scoring scheme
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -1

        # Perform the alignment

        alignments = aligner.align(seq1, seq2)

        # Get the best alignment
        best_alignment = alignments[0]

        # Extract the individual text from the best alignment
        alignment_str = str(best_alignment)
        alignment_lines = alignment_str.split('\n')
        target_seq = ''
        middle_line = ''
        query_seq = ''
        for i in range(0, len(alignment_lines), 4):
            target_line = alignment_lines[i]
            middle_line_parts = alignment_lines[i + 1].split()
            query_line = alignment_lines[i + 2]

            # Remove extraneous text from target and query lines
            target_line_parts = target_line.split()
            target_seq += target_line_parts[-1]
            target_seq = target_seq.replace("-", "X")  # reference
            query_line_parts = query_line.split()
            query_seq += query_line_parts[-1]
            query_seq = query_seq.replace("-", "X")  # read

            # Remove extraneous text from middle line
            middle_line += ' ' * (len(target_line_parts[0]) - len(middle_line_parts[0]))
            middle_line += middle_line_parts[-1]

        # =========== AI code exiting ===========

        # finding areas of insertion
        indexes = [i for i, letter in enumerate(target_seq) if letter == "X"]
        insertions = [query_seq[i] for i in indexes]

        # remove insertions from original read
        read = "".join(multi_delete(list(query_seq), indexes))

        # remove any trailing Xs
        for i in range(len(read) - 1, 0, -1):
            if not (read[i] == 'X' or (read[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])):
                read = read[:i + 1]
                break

        ref = seq1[:len(read)]

        print(read)
        print(ref)
        # finding snp positions (i = local | i + start_pos = global)
        for i in range(len(read)):
            if ref[i] != read[i]:
                print(ref[i], read[i])
                snp_list.append(SNP(i + start_pos, i))
                if not (i + start_pos in snp_compare_list):
                    DNA.global_snp.append(SNP(i + start_pos))
                    snp_compare_list.append(i + start_pos)
        print(len(read))
        # changing to DNA class
        read = np.array(list(read))
        print(len(read))

        final_sam = DNA(read, snp_list, start_pos, insertions, indexes)
        mod_sam_list = np.append(mod_sam_list, final_sam)

        snp_list = []

    return mod_sam_list


def reference_reader(filename):
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    l = []
    for b in fasta_sequences:
        l.append(b)
    return l[0]


def sam_reader(filename, referencefilename):
    samfile = pysam.AlignmentFile(filename, 'r')
    reference = reference_reader(referencefilename)

    index_pos = []
    reads = []

    for sam in samfile.fetch():
        # not using reads with soft or hard clips cause f*** you
        if not sam.cigartuples[0][0] > 2 or sam.cigartuples[-1][0] > 2:
            index_pos.append(sam.pos)
            reads.append(sam.seq)
        elif sam.cigartuples[0][0] > 2 and not sam.cigartuples[-1][0] > 2:
            index_pos.append(sam.pos)
            read = sam.seq
            read = read[sam.cigartuples[0][1]:]
            reads.append(read)
        elif not sam.cigartuples[0][0] > 2 and sam.cigartuples[-1][0] > 2:
            index_pos.append(sam.pos)
            read = sam.seq
            read = read[:sam.cigartuples[-1][1] * -1]
            reads.append(read)
        else:
            index_pos.append(sam.pos)
            read = sam.seq
            read = read[sam.cigartuples[0][1]:sam.cigartuples[-1][1] * -1]
            reads.append(read)

    new_reads = sam_aligner(reference.seq, reads, index_pos)

    return new_reads


def main():
    new_reads = sam_reader("SamFile.sam", "ReferenceFile.fasta")
    DNA.snp_purger(new_reads)
    log("debug", DNA.global_snp)


if __name__ == '__main__':
    main()
