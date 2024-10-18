#!python
# Annotate a set of J sequences, finding the frame alignment and index of the first nucleotide of the conserved PHE or TRP

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
import json

from receptor_utils import aux_formats


def main():
    args = get_parser().parse_args()

    if args.seq_file.lower().endswith('.json'):
        with open(args.seq_file, 'r') as fo:
            ref = json.load(fo)
        aux_formats.aux_from_json(ref, args.out_file, args.verbose)
    elif args.seq_file.lower().endswith('.fasta'):
        aux_formats.aux_from_fasta(args.seq_file, args.out_file, args.verbose)
    else:
        print("Error: reference file must be a .json or .fasta file")
        exit(1)


def get_parser():
    parser = argparse.ArgumentParser(
        description='Annotate a set of J sequences, finding the frame alignment and 1-based index of the first nucleotide of the conserved PHE or TRP')
    parser.add_argument('seq_file', help='J sequences to annotate (fasta or AIRR-C JSON format)')
    parser.add_argument('out_file', help='Annotation output in aux format for use with IgBlast')
    parser.add_argument('--verbose', '-v', help='Verbose output', action='store_true')
    return parser




if __name__ == "__main__":
    main()
