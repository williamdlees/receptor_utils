#!python
# Make igblast ndm file from a gapped reference set

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
import json
from receptor_utils import aux_formats


def get_parser():
    parser = argparse.ArgumentParser(description='Make IMGT-delineated igblast ndm file from an IMGT-gapped variable-region reference set or AIRR-C JSON format file.')
    parser.add_argument('ref_file', help='gapped reference set (.fasta, .json)')
    parser.add_argument('chain', help='chain as required by igblast (eg VH)')
    parser.add_argument('ndm_file', help='output file (ndm)')
    parser.add_argument('--cdr_coords', '-c', help='comma-separated list of 1-based co-ordinates for CDR1, CDR2 and CDR3-start (default 79,114,166,195,313)')
    return parser


def main():
    args = get_parser().parse_args()

    allowed_chains = ['VH', 'VK', 'VL', 'VA', 'VB', 'VD', 'VG']

    if args.chain not in allowed_chains:
        print(f"Error: chain must be one of {', '.join(allowed_chains)}")
        exit(1)

    cdr_coords = [79, 114, 166, 195, 313]
    if args.cdr_coords is not None:
        if not args.ref_file.lower().endswith('.fasta'):
            print("Error: CDR coordinates are only used with fasta files")
            exit(1)
        cdr_coords = args.cdr_coords.split(',')
        if len(cdr_coords) != 5:
            print("Error: CDR coordinates must be a comma-separated list of 5 integers")
            exit(1)
        try:
            cdr_coords = [int(x) for x in cdr_coords]
        except:
            print("Error: CDR coordinates must be a comma-separated list of 5 integers")
            exit(1)

        last_coord = 0
        for coord in cdr_coords:
            if coord < 0 or coord <= last_coord:
                print("Error: CDR coordinates must be positive integers in ascending order")
                exit(1)
            last_coord = coord

    if args.ref_file.lower().endswith('.json'):
        with open(args.ref_file, 'r') as fo:
            ref = json.load(fo)

        aux_formats.ndm_from_json(ref, args.ndm_file, args.chain, 'IMGT')
    elif args.ref_file.lower().endswith('.fasta'):
        aux_formats.ndm_from_fasta(args.ref_file, args.ndm_file, args.chain, cdr_coords)
    else:
        print("Error: reference file must be a .json or .fasta file")
        exit(1)


if __name__ == '__main__':
    main()
