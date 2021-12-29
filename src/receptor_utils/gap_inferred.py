# Gap sequences inferred by IgDiscover

import number_ighv
import argparse
import simple_bio_seq as simple

parser = argparse.ArgumentParser(description='Gap inferred sequences')
parser.add_argument('inferred_file', help='ungapped inferred sequences (fasta)')
parser.add_argument('ref_file', help='gapped reference set (fasta)')
parser.add_argument('out_file', help='output file containing gapped inferred sequences')
args = parser.parse_args()

inferred = simple.read_fasta(args.inferred_file)
refs = simple.read_fasta(args.ref_file)
ungapped_refs = {}

for id, seq in refs.items():
    ungapped_refs[id] = seq.replace('.', '')

gapped = {}

for id, seq in inferred.items():
    if '.' not in seq:
        gapped[id], aa, notes = number_ighv.gap_sequence(seq, refs, ungapped_refs)
        if notes is not None and len(notes) > 0:
            print('%s: %s' % (id, notes))
    else:
        gapped[id] = seq

simple.write_fasta(gapped, args.out_file)

