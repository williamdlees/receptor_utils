#!python
# Extract reference files for nominated species
# The current IMGT reference file can be downloaded from http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
from receptor_utils import simple_bio_seq as simple


def main():
    parser = argparse.ArgumentParser(description='Extract reference files for nominated species')
    parser.add_argument('imgt_file', help='gapped imgt reference file')
    parser.add_argument('species_name', help='species name for IMGT file in quotes, e.g. "Homo sapiens", "Macaca mulatta"')
    args = parser.parse_args()

    segs = ['IGHV', 'IGHD', 'IGHJ', 'CH']
    refs = simple.read_imgt_fasta(args.imgt_file, [args.species_name], segs)

    simple.write_fasta(refs[args.species_name]['IGHV'], '%s_IGHV_gapped.fasta' %args.species_name.replace(' ', '_'))

    ungapped = {}

    for seg in segs:
        ungapped[seg] = {}
        for id, seq in refs[args.species_name][seg].items():
            ungapped[seg][id] = seq.replace('.', '')

        simple.write_fasta(ungapped[seg], '%s_%s.fasta' % (args.species_name.replace(' ', '_'), seg))



if __name__ == "__main__":
    main()
