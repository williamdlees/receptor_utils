# Check the aux files against igblast

import argparse
from receptor_utils import simple_bio_seq as simple

parser = argparse.ArgumentParser(
    description='Check the aux files against igblast')
parser.add_argument('aux_file', help='Aux file to check')
parser.add_argument('igblast_file', help='Igblast file to check against')
args = parser.parse_args()


def read_whitespace_with_defaults(file_path, keys=None):
    records = []
    with open(file_path, 'r') as file:
        for line in file:
            elements = line.strip().split()  # Split by any whitespace
            if keys:
                elements += [''] * (len(keys) - len(elements))  # Fill missing elements with empty strings
                record = dict(zip(keys, elements))
            else:
                record = {f'col{i}': elem for i, elem in enumerate(elements)}
            records.append(record)
    return records


# Define the keys based on the expected columns
keys = ['#name', 'j_codon_frame', 'chain_type', 'j_cdr3_end', 'extra_bps']

aux_recs = simple.read_csv(args.aux_file, delimiter='\t')
igblast_recs = read_whitespace_with_defaults(args.igblast_file, keys=keys)

igblast = {x['#name']: x for x in igblast_recs}

for rec in aux_recs:
    if rec['#name'][0] == '#':
        continue

    if '*' not in rec['#name']:
        rec['#name'] = rec['#name'].split('-')[0] + '0-' + rec['#name'].split('-')[1] + '*00'

    if rec['#name'] not in igblast:
        print(f"Error: {rec['#name']} not found in igblast")
        continue

    if rec['j_codon_frame'] != igblast[rec['#name']]['j_codon_frame']:
        print(f"Error: {rec['#name']} j_codon_frame {rec['j_codon_frame']} != {igblast[rec['#name']]['j_codon_frame']}")
    if rec['chain_type'] != igblast[rec['#name']]['chain_type']:
        print(f"Error: {rec['#name']} chain_type {rec['chain_type']} != {igblast[rec['#name']]['chain_type']}")
    if rec['j_cdr3_end'] != igblast[rec['#name']]['j_cdr3_end']:
        print(f"Error: {rec['#name']} j_cdr3_end {rec['j_cdr3_end']} != {igblast[rec['#name']]['j_cdr3_end']}")
    if rec['extra_bps'] != igblast[rec['#name']]['extra_bps']:
        print(f"Error: {rec['#name']} extra_bps {rec['extra_bps']} != {igblast[rec['#name']]['extra_bps']}")