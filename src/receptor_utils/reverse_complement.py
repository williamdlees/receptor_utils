# Copyright (c) 2022 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


from receptor_utils import simple_bio_seq as simple
import argparse


def main():
    parser = argparse.ArgumentParser(description='Return the reverse complement of a nucleotide sequence')
    parser.add_argument('sequence', help='input sequence')
    args = parser.parse_args()

    rc = simple.reverse_complement(args.sequence)
    print(rc)
