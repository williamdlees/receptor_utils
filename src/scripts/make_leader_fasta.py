# Create a fasta file containing leader sequences for all V alleles in a given JSON file.

import argparse
import json
from receptor_utils import aux_formats as aux


def get_parser():
    parser = argparse.ArgumentParser(description="Create a fasta file containing leader sequences for all V alleles in a given JSON file.")
    parser.add_argument("json_file", help="Path to the input JSON file containing germline set information.")
    parser.add_argument("out_file", help="Path to the output FASTA file to be created.")
    parser.add_argument("--aa", action="store_true", help="Translate leader sequences to amino acids.")
    return parser


if __name__ == "__main__":
    
    args = get_parser().parse_args()
    
    with open(args.json_file, 'r') as f:
        germline_set = json.load(f)
    
    success = aux.write_leaders_to_fasta(germline_set, args.out_file, aa=args.aa)
    
    if success:
        print(f"Leader FASTA file successfully created at {args.out_file}.")
    else:
        print("No leader sequences were found. No FASTA file was created.")
