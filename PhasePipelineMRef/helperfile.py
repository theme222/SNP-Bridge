"""
TechnicalTools.py + SamReader.py
"""

from inspect import currentframe, getframeinfo
import os
import pysam
from Bio import SeqIO
from random import randint, choice
from Bio.Align import PairwiseAligner
import numpy as np
from collections import Counter


# Did you know you could do this? : 😄


def log(logtype, message, var=None):
    """
    like the logging import but less ugly
    :param var:
    :param message:
    :param logtype:
    """
    linenumber = currentframe().f_back.f_lineno
    if os.name == 'nt':
        filename = getframeinfo(currentframe().f_back).filename.split('\\')[-1]
    else:
        filename = getframeinfo(currentframe().f_back).filename.split('/')[-1]
    if logtype == 'info':
        print_message = f'\033[1;34m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    elif logtype == 'debug':
        print_message = f'\033[1;35m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    elif logtype == 'error':
        print_message = f'\033[0;31m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    elif logtype == 'critical':
        print_message = f'\033[1;31m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    else:
        print_message = f'\033[32m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    if var is None:
        print(print_message[:-10])
    else:
        print(print_message)


def time_convert(seconds):
    """
    turns absurd amounts of seconds into understandable text
    :param seconds:
    :return:
    """

    year = 31557600
    month = 2629800
    week = 604800
    day = 86400
    hour = 3600
    minute = 60
    second = 1
    millisecond = 0.001
    microsecond = 0.000001
    nanosecond = 0.000000001

    return_text = ''

    amount = seconds // year
    if amount > 0:
        return_text = return_text + f' {int(amount)} Year'
        seconds -= year * amount

    amount = seconds // month
    if amount > 0:
        return_text = return_text + f' {int(amount)} Month'
        seconds -= month * amount

    amount = seconds // week
    if amount > 0:
        return_text = return_text + f' {int(amount)} Week'
        seconds -= week * amount

    amount = seconds // day
    if amount > 0:
        return_text = return_text + f' {int(amount)} Day'
        seconds -= day * amount

    amount = seconds // hour
    if amount > 0:
        return_text = return_text + f' {int(amount)} Hour'
        seconds -= hour * amount

    amount = seconds // minute
    if amount > 0:
        return_text = return_text + f' {int(amount)} Minute'
        seconds -= minute * amount

    amount = seconds // second
    if amount > 0:
        return_text = return_text + f' {int(amount)} Second'
        seconds -= second * amount

    amount = int(seconds / millisecond)
    if amount > 0:
        return_text = return_text + f' {amount} MilliSecond'
        seconds -= millisecond * amount

    amount = int(seconds / microsecond)
    if amount > 0:
        return_text = return_text + f' {amount} MicroSecond'
        seconds -= microsecond * amount

    amount = int(seconds / nanosecond)
    if amount > 0:
        return_text = return_text + f' {amount} NanoSecond'
        seconds -= nanosecond * amount

    return return_text + ' '


def color_gen(style='rgb'):
    """ Makes a random color depending on style """
    color_names = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'cyan', 'black',
                   'magenta', 'yellow', 'aqua', 'purple', 'pink']
    if style == 'rgb' or style == 'RGB':
        return [randint(0, 255), randint(0, 255), randint(0, 255)]
    elif style == 'hex':
        return '#' + format(randint(0, 16777215), 'x')
    elif style == 'name' or style == 'color':
        return choice(color_names)


class TextFile:
    """Simple class to format text files"""

    def __init__(self):
        pass

    def __repr__(self):
        return 'Why would you do this?'

    @staticmethod
    def delnewline(filename, char_amount):
        """Deletes a set number of characters after a new line"""
        finishedtext = ''
        for line in open(filename, 'r'):
            finishedtext += line[char_amount:]
        filename_new = filename.split('.')[0] + '_new.txt'
        try:
            open(filename_new, 'x').write(finishedtext)
        except FileExistsError:
            log('error', 'Created file already exists', filename_new)

    @staticmethod
    def replaceword(filename, word, replacement=''):
        """finds the word and replaces with another thing"""
        finishedtext = ''
        for line in open(filename, 'r'):
            finishedtext += line.replace(word, replacement)
        filename_new = filename.split('.')[0] + '_new.txt'
        try:
            open(filename_new, 'x').write(finishedtext)
        except FileExistsError:
            log('error', 'Created file already exists', filename_new)


def multi_delete(list_, args):
    args = set(args)
    indexes = sorted(args, reverse=True)
    for index in indexes:
        del list_[index]
    return list_


class SNP:
    """ class to store SNP data type containing allele and global position"""

    def __init__(self, global_start_pos, local_start_pos=None, values=None):
        self.allele = []
        self.global_start_pos = global_start_pos
        if local_start_pos is None:  # used in read.snp
            self.local_start_pos = global_start_pos
        else:
            self.local_start_pos = local_start_pos
        if values is None:  # used in global snp
            self.values = []
        else:
            self.values = values

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

    def __init__(self, read=np.array([], dtype=str), snp=None, start_pos=0, insertions=None):
        if snp is None: snp = []
        if insertions is None:  insertions = np.array([], dtype=object)
        self.read = read  # Numpy array
        self.snp = snp  # Normal Python list containing SNP objects?
        self.start_pos = start_pos  # Int
        self.insertid = 0  # Int
        self.insertions = insertions  # Numpy array of insertion object

    def __repr__(self):
        return f'"{self.start_pos}"'
        #  return f'{list(self.read)} \n mutates at {self.snp}'

    """
    Adenine == 1
    Thymine == 2
    Cytosine == 3
    Guanine == 4
    """

    @property
    def size(self):
        return len(self.read)

    @property
    def read_with_ins(self):
        if len(self.insertions) != 0:
            read_with_in = self.read.copy()
            _temp = 0
            for index, ins in enumerate(self.insertions):
                read_with_in = np.insert(read_with_in, ins.local_start_pos, ins.values)
            return "".join(read_with_in)
        return "".join(self.read)

    @classmethod
    def snp_purger(cls, read_list):  # finding snps that aren't snps

        #  Purge useless snps
        snp_values = []
        final_snp_list = []
        for index, snp in enumerate(cls.global_snp):
            for index1, read in enumerate(read_list):
                if read.start_pos <= snp.global_start_pos <= read.start_pos + read.size - 1:  # if snp is in read
                    # Append the snp to the reads that didn't get detected as snps because it was teh same as the ref
                    if not snp.snp_checker(read.snp)[0]:
                        _temp = SNP(snp.global_start_pos)
                        _temp.local_start_pos = _temp.global_start_pos - read.start_pos
                        read.snp.append(_temp)
                        read_list[index1] = read
                    snp_values.append(read.read[snp.global_start_pos - read.start_pos])

            # Combine all words together
            chars = "".join(snp_values)
            # Count all characters
            counter = Counter(chars)
            # Get frequency and put in list
            total_occ = [(char, occ) for char, occ in counter.most_common()]
            if len(total_occ) >= 2:  # only include a snp when it has two values or more
                if total_occ[1][1] / total_occ[0][1] > 0.4:  # only if the two values are around the same amount
                    snp.values = [total_occ[0][0], total_occ[1][0]]  # gives the two possible values that the snp can be
                    final_snp_list.append(snp)

            snp_values = []

        cls.global_snp = final_snp_list
        cls.global_snp.sort(key=SNP.key)
        log('critical', cls.global_snp)

    @staticmethod
    def key(read):
        return read.start_pos


class Insertion(SNP):
    def __init__(self, global_insertion_pos=0, local_insert_pos=0, insertion=""):
        super().__init__(global_insertion_pos, local_insert_pos, values=insertion)

    '''
    AGTTAGGACAIIIATGATAGCCCAGTAGCA
              ^^^
    0123456789
    local = 10
    s[:10] = AGTTAGGACA
    '''


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

        # making reference that is slightly bigger that target sequence, so it can account for ins/del
        seq1 = reference[start_pos:start_pos + len(sam_seq) + 50]

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

        read = query_seq
        # remove any trailing Xs
        for i in range(len(read) - 1, 0, -1):
            if not (read[i] == 'X' or (read[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])):
                read = read[:i + 1]
                break

        ref = seq1[:len(read)]
        # finding areas of insertion
        indexes = [i for i, letter in enumerate(ref) if letter == "X"]
        insertions = [Insertion(ins + start_pos, ins, sam_seq[ins]) for ins in indexes]  # Make list of insertions

        # remove insertions from original read
        read = "".join(multi_delete(list(read), indexes))
        print(len(ref),len(read))
        # finding snp positions (i = local | i + start_pos = global)
        if len(ref) == len(read):
            for i in range(len(read)):
                if ref[i] != read[i]:
                    #  print(ref[i], read[i])
                    snp_list.append(SNP(i + start_pos, i))
                    if not (i + start_pos in snp_compare_list):
                        DNA.global_snp.append(SNP(i + start_pos))
                        snp_compare_list.append(i + start_pos)

            # changing to DNA class
            read = np.array(list(read))

            final_sam = DNA(read, snp_list, start_pos, np.array(insertions))
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

    for sam in (samfile.fetch()):
        # cutting the reads that have soft clips and hard clips out
        print(sam.cigartuples)
        if not sam.cigartuples[0][0] > 2 or sam.cigartuples[-1][0] > 2:
            index_pos.append(sam.pos)
            reads.append(sam.seq)
        elif sam.cigartuples[0][0] > 2 and not sam.cigartuples[-1][0] > 2:
            index_pos.append(sam.pos + sam.cigartuples[0][1])
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

    return new_reads, reference.seq


def remove_unaligned_reads(input_file, output_file):
    with pysam.AlignmentFile(input_file, "r") as input_sam, \
            pysam.AlignmentFile(output_file, "w", template=input_sam) as output_sam:
        for read in input_sam:
            if read.cigartuples is not None:
                output_sam.write(read)
    os.remove(input_file)