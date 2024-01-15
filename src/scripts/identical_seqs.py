#!python
# List identical sequences and sub-sequences in a fasta file

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
from receptor_utils import simple_bio_seq as simple


def get_parser():
    parser = argparse.ArgumentParser(description='List identical sequences and sub-sequences in a fasta file')
    parser.add_argument('input_file', help='gapped imgt reference file')
    return parser


def main():
    args = get_parser().parse_args()

    seqs = simple.read_fasta(args.input_file)

    seq_index = {}

    for seq_name, seq in seqs.items():
        if seq not in seq_index:
            seq_index[seq] = []
        seq_index[seq].append(seq_name)

    print('Identical sequences:')

    for names in seq_index.values():
        if len(names) > 1:
            print(', '.join(names))

    print('Sub-sequences:')
    for n1, s1 in seqs.items():
        for n2, s2 in seqs.items():
            if s1 != s2  and s1 in s2:
                print(f'{n1} ({s1})is a sub-sequence of {n2} ({s2})')



if __name__ == "__main__":
    main()
