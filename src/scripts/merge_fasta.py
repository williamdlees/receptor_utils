# Merge two fasta files, keeping records in file 2 that do not match sequences in file 1

import argparse
from receptor_utils import simple_bio_seq as simple


def get_parser():
    parser = argparse.ArgumentParser(description='Merge two fasta files, keeping records in file 2 that do not match sequences in file 1')
    parser.add_argument('file_1', help='first file to merge')
    parser.add_argument('file_2', help='second file to merge')
    parser.add_argument('outfile', help='output file')
    return parser


def main():
    args = get_parser().parse_args()

    seqs = simple.read_fasta(args.file_1)
    seqs_map = {v: k for k, v in seqs.items()}

    merge_seqs = simple.read_fasta(args.file_2)

    added = 0

    for m_name, m_seq in merge_seqs.items():
        if m_seq not in seqs_map:
            seqs[m_name] = m_seq
            added += 1

    simple.write_fasta(seqs, args.outfile)
    print('%d items merged from %s.' % (added, args.file_2))


if __name__ == "__main__":
    main()
