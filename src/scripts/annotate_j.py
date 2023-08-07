#!python
# Annotate a set of J sequences, finding the frame alignment and index of the first nucleotide of the conserved PHE or TRP

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
from receptor_utils import simple_bio_seq as simple
import argparse
import re
from collections import namedtuple
import json


def main():
    args = get_parser().parse_args()

    if args.seq_file.lower().endswith('.json'):
        seqs = seqs_from_json(args)
    elif args.seq_file.lower().endswith('.fasta'):
        seqs = simple.read_fasta(args.seq_file)
    else:
        print("Error: reference file must be a .json or .fasta file")
        exit(1)

    
    annotations = []
    
    for seq_name, seq in seqs.items():
        solutions = []
        Result = namedtuple('Result', ['nt_seq', 'frame', 'cdr3_pos', 'trans', 'extra_bps'])
        for frame in range(3):
            trans = simple.translate(seq[frame:])

            if '*' in trans:
                continue

            for m in re.finditer('[WF]G.G', trans):
                solutions.append(Result._make([
                    seq,
                    frame,
                    frame + 3*m.start(),
                    trans,
                    (len(seq) - frame) % 3
                ]))

        if solutions:
            if len(solutions) > 1:
                print('WARNING: multiple solutions for {seq_name}')
                
            for solution in solutions:
                print(f'{seq_name}: j_codon_frame {solution.frame} j_cdr3_end {solution.cdr3_pos-1} exta bps {solution.extra_bps}')

                loc_count = '      '
                for num in range(len(solution.trans)):
                    loc_count += f' {num+1:<2}'
                print(loc_count)

                aa = '      '
                for num in range(len(solution.trans)):
                    aa += f' {solution.trans[num]} '
                print(aa)

                nt = f'  {solution.nt_seq[:solution.frame]:>3} {solution.nt_seq[solution.frame:]}'
                print(nt + '\n\n')

                annotations.append({
                    '#name': seq_name,
                    'j_codon_frame': solution.frame,
                    'chain_type': f'{seq_name[3]}{seq_name[2]}',
                    'j_cdr3_end': solution.cdr3_pos-1,
                    'extra_bps': solution.extra_bps,
                })
        else:
            print(f'{seq_name}: no solutions')

    simple.write_csv(args.out_file, annotations, delimiter='\t')


def get_parser():
    parser = argparse.ArgumentParser(
        description='Annotate a set of J sequences, finding the frame alignment and 1-based index of the first nucleotide of the conserved PHE or TRP')
    parser.add_argument('seq_file', help='J sequences to annotate (fasta or AIRR-C JSON format)')
    parser.add_argument('out_file', help='Annotation output in aux format for use with IgBlast')
    return parser


def seqs_from_json(args):
    with open(args.seq_file, 'r') as fo:
        ref = json.load(fo)

    if isinstance(ref["GermlineSet"], list):
        germline_set = ref["GermlineSet"][0]
    else:
        germline_set = ref["GermlineSet"]

    seqs = {}

    for allele in germline_set['allele_descriptions']:
        if allele['sequence_type'] == 'J':
            seqs[allele['label']] = allele['coding_sequence']

    return seqs


if __name__ == "__main__":
    main()
