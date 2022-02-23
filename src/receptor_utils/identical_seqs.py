#!python
# List identical sequences and sub-sequences in a fasta file

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
from receptor_utils import simple_bio_seq as simple


def main():
    parser = argparse.ArgumentParser(description='List identical sequences and sub-sequences in a fasta file')
    parser.add_argument('input_file', help='gapped imgt reference file')
    args = parser.parse_args()

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
    for n1 in seqs:
        for n2 in seqs:
            if n1 != n2 and seqs[n1] != seqs[n2] and n1 in n2:
                print('%s is a sub-sequence of %s' % (n1, n2))



if __name__ == "__main__":
    main()
