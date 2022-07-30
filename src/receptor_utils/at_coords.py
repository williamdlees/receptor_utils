#!python
# Extract the sequence at the specified 1-based co-ordinates in the reference sequence

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
from receptor_utils import simple_bio_seq as simple


def main():
    parser = argparse.ArgumentParser(description=' Extract the sequence at the specified 1-based co-ordinates in the file')
    parser.add_argument('input_file', help='File containing the reference sequence (FASTA)')
    parser.add_argument('start', help='Start co-ord (1-based)', type=int)
    parser.add_argument('end', help='End co-ord (1-based)', type=int)
    parser.add_argument('-r', help='Reverse-complement the result', action='store_true')
    args = parser.parse_args()

    if args.end < args.start or args.end < 1 or args.start < 1:
        quit()

    ref = simple.read_single_fasta(args.input_file)
    seq = ref[args.start-1:args.end]

    if args.r:
        seq = simple.reverse_complement(seq)

    print(seq)




