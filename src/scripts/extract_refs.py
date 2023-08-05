#!python
# Extract reference files for nominated species
# The current IMGT reference file can be downloaded from https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
from receptor_utils import simple_bio_seq as simple
import urllib.request
import io


def get_parser():
    parser = argparse.ArgumentParser(description='Extract reference files for nominated species from a gapped IMGT reference file')
    parser.add_argument('species_name', help='species name for IMGT file in quotes, e.g. "Homo sapiens", "Macaca mulatta"')
    parser.add_argument('-r', '--ref', help='gapped imgt reference file. If not specified, the current file will be downloaded from http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP')
    parser.add_argument('-L', '--locus', help='locus identifier, e.g. "IGH" (default), "IGK", "TRB", "TRG"', default='IGH', required=False)
    parser.add_argument('-F', '--functional_only', action='store_true')
    return parser

def main():
    args = get_parser().parse_args()

    imgt_url = "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP"

    segs = {
        'IGH': ['IGHV', 'IGHD', 'IGHJ', 'CH'],
        'IGK': ['IGKV', 'IGKJ', 'IGKC'],
        'IGL': ['IGLV', 'IGLJ', 'IGLC'],
        'TRA': ['TRAV', 'TRAJ', 'TRAC'],
        'TRB': ['TRBV', 'TRBD', 'TRBJ', 'TRBC'],
        'TRG': ['TRGV', 'TRGJ', 'TRGC'],
        'TRD': ['TRDV', 'TRDD', 'TRDJ', 'TRDC']
    }

    if args.ref:
        refs = simple.read_imgt_fasta(args.ref, [args.species_name], segs[args.locus], functional_only=args.functional_only)
    else:
        with urllib.request.urlopen(imgt_url) as fi:
            refs = simple.read_imgt_fasta(io.StringIO(fi.read().decode('utf-8')), [args.species_name], segs[args.locus], functional_only=args.functional_only)

    simple.write_fasta('%s_IGHV_gapped.fasta'.replace('IGHV', args.locus + 'V') %args.species_name.replace(' ', '_'), refs[args.species_name][args.locus + 'V'])

    ungapped = {}

    for seg in segs[args.locus]:
        ungapped[seg] = {}
        for id, seq in refs[args.species_name][seg].items():
            ungapped[seg][id] = seq.replace('.', '')

        simple.write_fasta('%s_%s.fasta' % (args.species_name.replace(' ', '_'), seg), ungapped[seg])


if __name__ == "__main__":
    main()
