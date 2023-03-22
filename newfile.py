from Bio.Align import PairwiseAligner
from TechnicalTools import multi_delete


class SamRead():
    def __init__(self, read='', start_pos=0, insertions=None, insertion_pos=None):
        if insertion_pos is None:  insertion_pos = []
        if insertions is None:  insertions = []

        self.read = read  # without insertions
        self.start_pos = start_pos
        self.insertions = insertions
        self.insertion_pos = insertion_pos

    def __repr__(self):
        return f"read: {self.read} insertions: {self.insertions} with position: {self.insertion_pos}"


def samAligner(reference, sam_list, sam_start_pos):
    """
    reference: String of reference gene
    sam_list: list of strings
    sam_start_pos: list of integers of starting position 0
    return = positions of snps , modified sam_list
    """

    # return variables
    snp_pos = []
    mod_sam_list = []

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
        for i in range(len(read), 0, -1):
            if read[i] != 'X':
                read = read[:i]
                break

        ref = seq2[:len(read)-1]

        # finding snp positions
        for i in range(len(read)):
            if ref[i] != read[i]:
                snp_pos.append(i)


        # changing to SamRead class
        final_sam = SamRead(read, start_pos, insertions, indexes)
        mod_sam_list.append(final_sam)

    return snp_pos,mod_sam_list
