#!python
# Make a mixcr-compatible json file from MIARR json file

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import argparse
import json
from receptor_utils import aux_formats


def get_parser():
    parser = argparse.ArgumentParser(description='Make a mixcr-compatible json file from MIARR json file.')
    parser.add_argument('miairr_ref_file', help='MIAIRR reference set (.json)')
    parser.add_argument('mixcr_ref_file', help='output file (.json)')
    parser.add_argument('--taxon_id', "-t", help='taxon ID for the reference set (e.g. 9606 for human)')
    parser.add_argument('--species_name', "-s", help='species name for the reference set (e.g. "Homo sapiens")')
    return parser


def main():
    args = get_parser().parse_args()

    with open(args.miairr_ref_file, 'r') as fo:
        ref = json.load(fo)

    if isinstance(ref["GermlineSet"], list):
        germline_set = ref["GermlineSet"][0]
    else:
        germline_set = ref["GermlineSet"]

    taxon_id = args.taxon_id

    if not taxon_id:
        if 'species' in germline_set and 'id' in germline_set['species']:
            taxon_id = germline_set['species']['id']
            if ':' in taxon_id:
                taxon_id = taxon_id.split(':')[-1]
        else:
            print("Taxon ID not found in the reference file - please specify on the command line.")
            return
        
    species_name = args.species_name

    if not species_name:
        if 'species' in germline_set and 'label' in germline_set['species']:
            species_name = germline_set['species']['label']
        else:
            print("Species name not found in the reference file - please specify on the command line.")
            return
        
    uri_prefix = f"file://{args.miairr_ref_file}"

    mixcr_ref = aux_formats.make_mixcr_ref(ref, uri_prefix, taxon_id, species_name)

    with open(args.mixcr_ref_file, 'w') as fo:
        json.dump(mixcr_ref, fo, indent=4)

if __name__ == '__main__':
    main()
