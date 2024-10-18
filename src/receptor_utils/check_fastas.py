# Check the aux files against igblast

import argparse
from receptor_utils import simple_bio_seq as simple

parser = argparse.ArgumentParser(
    description='Check alignment of sequences in two fasta files')
parser.add_argument('file1', help='file to check')
parser.add_argument('file2', help='reference to check against')
args = parser.parse_args()

f1 = simple.read_fasta(args.file1)
f2 = simple.read_fasta(args.file2)

for name, seq1 in f1.items():
    if name not in f2:
        print(f"{name} not found in reference")
        continue

    seq2 = f2[name]
    if len(seq1) > len(seq2):
        print(f"{name} sequence longer than reference")
        seq1 = seq1[:len(seq2)]
    elif len(seq1) < len(seq2):
        print(f"{name} sequence shorter than reference")
        seq2 = seq2[:len(seq1)]

    if seq1 != seq2:
        print(f"{name} sequence mismatch")