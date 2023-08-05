#!python
# Gap sequences inferred by IgDiscover

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


from receptor_utils import number_v
import argparse
from receptor_utils import simple_bio_seq as simple
import Bio.Data.CodonTable


def get_parser():
    parser = argparse.ArgumentParser(description='Gap the sequences in a FASTA file, using the closest sequence in a gapped reference set as template')
    parser.add_argument('inferred_file', help='ungapped inferred sequences (fasta)')
    parser.add_argument('ref_file', help='gapped reference set (fasta)')
    parser.add_argument('out_file', help='output file containing gapped inferred sequences')
    return parser


def main():
    args = get_parser().parse_args()

    inferred = simple.read_fasta(args.inferred_file)
    refs = simple.read_fasta(args.ref_file)
    
    # Check that reference sequences have conserved residues - remove any that don't
    
    for ref, seq in list(refs.items()):
        try:
            # if the sequence has a leading gap, make it end on a codon boundary for the purpose of the checks

            trial_seq = seq.replace('.', '-')

            if trial_seq[0] == '-':
                aa = simple.translate(trial_seq, ignore_partial_codon=True)
            else:
                aa = simple.translate(trial_seq, ignore_partial_codon=False)

            notes = number_v.check_conserved_residues(aa)

            if notes:
                if 'truncated' in notes:
                    print(f"Warning: reference sequence {ref}: {notes} (if you don't want to use this as a reference sequence, please remove it from the file)")
                else:
                    print(f'Removing {ref} from reference: {notes}')
                    del refs[ref]
        except Bio.Data.CodonTable.TranslationError as e:
            print(f'Removing {ref} from reference: {e}\n{aa}')
            del refs[ref]

    ungapped_refs = {}

    for id, seq in refs.items():
        ungapped_refs[id] = seq.replace('.', '')

    gapped = {}

    for id, seq in inferred.items():
        gapped[id], aa, notes = number_v.gap_sequence(seq, refs, ungapped_refs)
        if notes:
            print('%s: %s' % (id, notes))

    simple.write_fasta(args.out_file, gapped)



if __name__ == "__main__":
    main()
