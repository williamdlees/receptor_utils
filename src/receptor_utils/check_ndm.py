# Check the ndm files against igblast

import argparse
from receptor_utils import simple_bio_seq as simple

parser = argparse.ArgumentParser(
    description='Check the ndm files against igblast')
parser.add_argument('ndm_file', help='ndm file to check')
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
keys = ['#name', 'FWR1 start', 'FWR1 stop', 'CDR1 start', 'CDR1 stop', 'FWR2 start', 'FWR2 stop', 'CDR2 start', 'CDR2 stop', 'FWR3 start', 'FWR3 stop', 'chain type', 'frame start']

ndm_recs = simple.read_csv(args.ndm_file, delimiter='\t')
igblast_recs = read_whitespace_with_defaults(args.igblast_file, keys=keys)

igblast = {x['#name']: x for x in igblast_recs}
error_count = 0

for rec in ndm_recs:
    if rec['#name'][0] == '#':
        continue

    if '*' not in rec['#name']:
        rec['#name'] = rec['#name'].split('-')[0] + '0-' + rec['#name'].split('-')[1] + '*00'

    igblast_key = rec['#name']
    if igblast_key not in igblast:
        igblast_key = igblast_key.split('*')[0] + '*01'
        if igblast_key not in igblast:
            print(f"{igblast_key} not found in igblast. *01 allele also not found: skipping")
            continue
        print(f"{rec['#name']} not found in igblast: using {igblast_key}")
        
    for coord in keys[1:]:
        if rec[coord] != igblast[igblast_key][coord]:
            print(f"Error: {rec['#name']} {coord} {rec[coord]} != {igblast[igblast_key][coord]}")
            error_count += 1

if error_count == 0:
    print('No errors found')
else:
    print(f'{error_count} errors found')
