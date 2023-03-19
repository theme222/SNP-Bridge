import json
import numpy as np
from Bio import SeqIO
from random import choice, randint


def fastA_parser(filename, count):
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    lst = []
    for fasta in fasta_sequences:
        lst.append(str(fasta.seq))
    return lst[:count]


def fastA_blender(fasta_reads, avg_depth=10, read_size=150):
    return_list = np.array([],dtype="<U4")
    for i in range(round(avg_depth * len(fasta_reads[0]) / (read_size * 2))):  # the amount of reads
        chosen_read = choice(fasta_reads)
        start_pos = randint(0, len(chosen_read) - (read_size * 2))
        return_list = np.append(return_list, chosen_read[start_pos:start_pos + (read_size * 2)])

    return return_list


def fastQ_maker(filename, reads):
    with open(f'{filename}.fastq', 'w') as file:
        for index, read in enumerate(reads):
            file.write(f"@SimultaedRead-{index}" + "\n")
            file.write(read + "\n")
            file.write("+" + "\n")
            file.write("~" * len(read) + "\n")

        file.close()


def main():
    with open("config.json") as json_data:
        configs = json.load(json_data)
    fasta_reads = fastA_parser("FastABase.fasta", 2)
    fastQ_maker("output", fastA_blender(fasta_reads,configs["Average Depth"],configs["Read Size"]))


if __name__ == '__main__':
    main()
