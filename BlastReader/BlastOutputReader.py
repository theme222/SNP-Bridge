import os
import re
import pandas as pd
from Bio.Blast.NCBIWWW import qblast
from TechnicalTools import log, time_convert
from Bio import SeqIO
import time
from Bio.Blast import NCBIXML


def blast(seq):
    try:
        blast_results = qblast('blastn', 'nr', seq, megablast=True, hitlist_size=10)
        # Parse the result using NCBIXML.read()
        blast_record = NCBIXML.read(blast_results)
        titles = []
        strings = ['HLA-A*', 'HLA-B*', 'HLA-C*', 'HLA-DPB1*', 'HLA-DQB1*', 'HLA-DRB1*']
        # Access information from the BLAST record
        for alignment in blast_record.alignments:
            print(alignment.title.split())
            for tiny_string in alignment.title.split():
                for specstring in strings:
                    if specstring in tiny_string:
                        _temp = tiny_string
                        if ':' in _temp:
                            final = ''
                            final += _temp.split(':')[0] + ':' + _temp.split(':')[1][:2]
                            titles.append(final)
        return set(titles)
    except KeyboardInterrupt:
        log('debug', f'KeyboardInterrupt')


def read_csv(file_path):
    df = pd.read_csv(file_path)
    a = df['hla_adr_allele']
    b = df['hla_adr_drug_id']
    c = df['hla_adr_cohort_ethnicity_id']
    my_dict = {key: (val1, val2) for key, val1, val2 in zip(a, b, c)}
    return my_dict


def find_in_dict(queries, dict):
    for query in queries:
        print(query)
        try:
            print(
                f'In gene {query} you have a risk of have an ADR against {dict[query][0]} if you are {dict[query][1]}')
        except KeyError:
            print('key not found')


def reference_reader(filename):
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    l = []
    for b in fasta_sequences:
        l.append(b)
    return l[0]


def main():
    starttime = time.time()
    fasta_seq = reference_reader('ADRFastA.fasta').seq
    log('info', 'Blasting ...')
    infos = blast(fasta_seq)
    data = read_csv('HLA-ADR - hla_adr.csv')
    find_in_dict(infos,data)

    log('debug', f'Total runtime of {time_convert(time.time() - starttime)}')


if __name__ == '__main__':
    main()
