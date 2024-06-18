#!python
# Extract the sequence at the specified 1-based co-ordinates in the reference sequence

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
from receptor_utils import simple_bio_seq as simple


def get_parser():
    parser = argparse.ArgumentParser(description=' Extract the sequence at the specified 1-based co-ordinates in the file')
    parser.add_argument('input_file', help='File containing the reference sequence (FASTA)')
    parser.add_argument('start', help='Start coord (1-based)', type=int)
    parser.add_argument('end', help='End coord (1-based)', type=int)
    parser.add_argument('-n', help='Label of sequence in file (must specify if there is more than one sequence in the file)')
    parser.add_argument('-r', help='Reverse-complement the result', action='store_true')
    parser.add_argument('-xr', help='Extract the result from the reverse-complement of the sequence (coords refer to the reverse complement sequence)', action='store_true')
    return parser


def main():
    args = get_parser().parse_args()

    if args.end < args.start or args.end < 1 or args.start < 1:
        quit()

    if not args.n:
        ref = simple.read_single_fasta(args.input_file)
    else:
        refs = simple.read_fasta(args.input_file)
        if args.n not in refs:
            print(f'Sequence {args.n} not found in file')
            print(f'Sequences in file: {", ".join(refs.keys())}')
            quit()
        ref = refs[args.n]

    if args.xr:
        ref = simple.reverse_complement(ref)

    seq = ref[args.start-1:args.end]

    if args.r:
        seq = simple.reverse_complement(seq)
        
    print(seq)




