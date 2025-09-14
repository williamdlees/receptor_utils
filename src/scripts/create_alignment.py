#!python
# Create a formatted alignment display from gapped sequences

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE

import argparse
from receptor_utils import simple_bio_seq as simple
from receptor_utils import sequence_alignment


def get_parser():
    parser = argparse.ArgumentParser(description='Create a formatted alignment display from gapped sequences')
    parser.add_argument('input_file', help='FASTA file containing gap-aligned nucleotide sequences')
    parser.add_argument('sequence_type', choices=['V', 'D', 'J'], help='Type of sequence (V, D, or J) - currently only V is implemented')
    parser.add_argument('output_file', help='Output file for the formatted alignment')
    parser.add_argument('--codon_wrap', '-w', type=int, default=20, help='Number of codons per line before wrapping (default: 20)')
    parser.add_argument('--v_coords', '-c', help='Comma-separated list of 1-based codon coordinates for CDR1, CDR2 and CDR3-start (default: 27,38,56,65,105)')
    parser.add_argument('--filter', '-f', help='Only include sequences in the file that contain this substring (e.g. IGHV1-18)')
    return parser


def main():
    args = get_parser().parse_args()
    
    # Parse V coordinates if provided
    v_coords = None
    if args.v_coords:
        try:
            coords = [int(x.strip()) for x in args.v_coords.split(',')]
            if len(coords) != 5:
                print("Error: v_coords must contain exactly 5 integers")
                exit(1)
            
            # Validate coordinates are ascending
            for i in range(1, len(coords)):
                if coords[i] <= coords[i-1]:
                    print("Error: v_coords values must be in ascending order")
                    exit(1)
            
            v_coords = tuple(coords)
        except ValueError:
            print("Error: v_coords must be a comma-separated list of integers")
            exit(1)
    
    # Read input sequences
    try:
        sequences = simple.read_fasta(args.input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        exit(1)

    # Apply filter if specified
    if args.filter:
        sequences = {k: v for k, v in sequences.items() if args.filter in k}
    
    if not sequences:
        if args.filter:
            print(f"Error: No sequences matching '{args.filter}' found in input file")
        else:
            print("Error: No sequences found in input file")
        exit(1)
    
    # Create alignment
    try:
        alignment = sequence_alignment.create_alignment(
            sequences,
            sequence_type=args.sequence_type,
            codon_wrap=args.codon_wrap,
            v_coords=v_coords
        )
    except Exception as e:
        print(f"Error creating alignment: {e}")
        exit(1)
    
    # Write output
    try:
        with open(args.output_file, 'w') as f:
            f.write(alignment)
        print(f"Alignment written to {args.output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        exit(1)


if __name__ == '__main__':
    main()
