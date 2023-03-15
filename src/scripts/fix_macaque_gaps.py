# Remove unwanted inserted codons in IMGT alignment of macaque IG sequences

from receptor_utils import simple_bio_seq as simple
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description='Remove unwanted inserted codons in IMGT alignment of rhesus macaque IG sequences')
    parser.add_argument('infile', help='gapped file from IMGT')
    parser.add_argument('outfile', help='gapped file confirming to the customary alignment (without inserted codons)')
    parser.add_argument('chain', help='IGH, IGK or IGL')
    return parser


unwanted_codons = {
    'IGL': [21, 52, 53],
    'IGK': [21],
    'IGH': [16, 28],
}


def main():
    parser = get_parser()
    args = parser.parse_args()

    if args.chain not in ['IGH', 'IGK', 'IGL']:
        print("'chain' must be one of IGH, IGK, IGL")
        exit(0)

    seqs = simple.read_fasta(args.infile)
    unwanted_seqs = []

    for name, seq in seqs.items():
        for unw in sorted(unwanted_codons[args.chain], reverse=True):
            if seq[(unw-1)*3] == '.':
                seq = seq[:(unw-1)*3] + seq[3 + (unw-1)*3:]
            else:
                print(f"Unwanted codon at {unw} in sequence {name} is non-blank: dropping sequence")
                unwanted_seqs.append(name)
                break
            seqs[name] = seq

    for unw in unwanted_seqs:
        del seqs[unw]

    simple.write_fasta(seqs, args.outfile)


if __name__ == "__main__":
    main()
