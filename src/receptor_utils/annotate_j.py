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


def main():
    parser = argparse.ArgumentParser(description='Annotate a set of J sequences, finding the frame alignment and 1-based index of the first nucleotide of the conserved PHE or TRP')
    parser.add_argument('seq_file', help='J sequences to annotate (fasta)')
    parser.add_argument('out_file', help='Annotation output (csv)')
    args = parser.parse_args()
    
    seqs = simple.read_fasta(args.seq_file)
    annotations = []
    
    for seq_name, seq in seqs.items():
        solutions = []
        Result = namedtuple('Result', ['nt_seq', 'frame', 'cdr3_pos', 'trans'])
        for frame in range(3):
            trans = simple.translate(seq[frame:])

            if 'X' in trans:
                continue

            for m in re.finditer('[WF]G.G', trans):
                solutions.append(Result._make([
                    seq,
                    frame,
                    frame + 3*m.start(),
                    trans
                ]))

        if solutions:
            if len(solutions) > 1:
                print('WARNING: multiple solutions for {seq_name}')
                
            for solution in solutions:
                print(f'{seq_name}: j_codon_frame {solution.frame+1} j_cdr3_end {solution.cdr3_pos+1}')

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
                    'name': seq_name,
                    'nt_sequence': seq,
                    'j_codon_frame': solution.frame+1,
                    'j_cdr3_end': solution.cdr3_pos+1,
                })
        else:
            print(f'{seq_name}: no solutions')

    simple.write_csv(args.out_file, annotations)

if __name__ == "__main__":
    main()
