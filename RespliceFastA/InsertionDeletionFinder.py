from Bio import SeqIO
from Bio.Align import PairwiseAligner


def fastA_parser(filename, count):
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    lst = []
    for fasta in fasta_sequences:
        lst.append(fasta)
    return lst[:count]


def sequenceAligner(seq1, seq2):  # AI Code lol

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
    target_seq = target_seq.replace("-", "X")
    query_line_parts = query_line.split()
    query_seq += query_line_parts[-1]
    query_seq = query_seq.replace("-", "X")

    # Remove extraneous text from middle line
    middle_line += ' ' * (len(target_line_parts[0]) - len(middle_line_parts[0]))
    middle_line += middle_line_parts[-1]

    return [target_seq, middle_line, query_seq]


def main():
    sequences = fastA_parser("FastAFile.fasta", 2)
    seq1, seq2 = sequences[0].seq, sequences[1].seq
    sequence = sequenceAligner(seq1, seq2)
    # Print the extracted text without numbers
    SNP_indices = [pos for pos, char in enumerate(sequence[1]) if char == '.']
    indel_indices = [pos for pos, char in enumerate(sequence[1]) if char == '-']
    print('Target seq    :', sequence[0])
    print('Middle line  : ', sequence[1])
    print('Query sequence:', sequence[2])
    print(SNP_indices)
    print(indel_indices)


if __name__ == '__main__':
    main()
