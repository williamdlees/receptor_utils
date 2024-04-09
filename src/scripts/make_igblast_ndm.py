#!python
# Make igblast ndm file from a gapped reference set

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
from receptor_utils import simple_bio_seq as simple
import argparse
import json

def get_parser():
    parser = argparse.ArgumentParser(description='Make IMGT-delineated igblast ndm file from an IMGT-gapped variable-region reference set or AIRR-C JSON format file.')
    parser.add_argument('ref_file', help='gapped reference set (.fasta, .json)')
    parser.add_argument('chain', help='chain as required by igblast (eg VH)')
    parser.add_argument('ndm_file', help='output file (ndm)')
    parser.add_argument('--cdr_coords', '-c', help='comma-separated list of 1-based co-ordinates for CDR1, CDR2 and CDR3-start (default 79,114,166,195,313)')
    return parser


def main():
    args = get_parser().parse_args()

    allowed_chains = ['VH', 'VK', 'VL', 'VA', 'VB', 'VD', 'VG']

    if args.chain not in allowed_chains:
        print(f"Error: chain must be one of {', '.join(allowed_chains)}")
        exit(1)

    cdr_coords = [79, 114, 166, 195, 313]
    if args.cdr_coords is not None:
        if not args.ref_file.lower().endswith('.fasta'):
            print("Error: CDR coordinates are only used with fasta files")
            exit(1)
        cdr_coords = args.cdr_coords.split(',')
        if len(cdr_coords) != 5:
            print("Error: CDR coordinates must be a comma-separated list of 5 integers")
            exit(1)
        try:
            cdr_coords = [int(x) for x in cdr_coords]
        except:
            print("Error: CDR coordinates must be a comma-separated list of 5 integers")
            exit(1)
        
        last_coord = 0
        for coord in cdr_coords:
            if coord < 0 or coord <= last_coord:
                print("Error: CDR coordinates must be positive integers in ascending order")
                exit(1)
            last_coord = coord
    
    if args.ref_file.lower().endswith('.json'):
        from_json(args)
    elif args.ref_file.lower().endswith('.fasta'):
        from_fasta(args, cdr_coords)
    else:
        print("Error: reference file must be a .json or .fasta file")
        exit(1)


def from_fasta(args, cdr_coords):
    seqs = simple.read_fasta(args.ref_file)
    
    
    # desired headers:
    # gene/allele name, FWR1 start, FWR1 stop, CDR1 start, CDR1 stop, FWR2 start, FWR2 stop, CDR2 start, CDR2 stop, FWR3 start, FWR3 stop, chain type, coding frame start.
    # FWR/CDR positions are 1-based while the coding frame start positions are 0-based
    
    omissions = False

    with open(args.ndm_file, 'w') as fo:
        for label, seq in seqs.items():

            # If the sequence is truncated at the 5' end, don't include the sequence

            if seq[0] == '.':
                print(f"Omitting allele {label} as it is truncated at the 5' end")
                omissions = True
                continue

            fo.write(record_from_gapped_seq(label, seq, args.chain, cdr_coords))

    if omissions:
        print('Omissions will not affect results provided at least one non-truncated allele of the gene is present.')


def record_from_gapped_seq(label, seq, chain, cdr_coords):
    # Python-based ranges of IMGT elements
    # default: 79,114,166,195,313
    # <---------------------------------- FR1-IMGT -------------------------------->______________ CDR1-IMGT ___________<-------------------- FR2-IMGT ------------------->___________ CDR2-IMGT ________<----------------------------------------------------- FR3-IMGT ----------------------------------------------------> CDR3-IMGT
    imgt_fr1 = (0, cdr_coords[0] - 1)
    imgt_cdr1 = (cdr_coords[0] - 1, cdr_coords[1])
    imgt_fr2 = (cdr_coords[1], cdr_coords[2] - 1)
    imgt_cdr2 = (cdr_coords[2] - 1, cdr_coords[3])
    imgt_fr3 = (cdr_coords[3], cdr_coords[4] - 1)

    rec = []
    rec.append(label)
    pos = 1
    rec.append(str(pos))     # FWR start
    pos += len(seq[slice(*imgt_fr1)].replace('.', '')) - 1
    rec.append(str(pos))    # FWR1 stop
    pos += 1
    rec.append(str(pos))     # CDR1 start
    pos += len(seq[slice(*imgt_cdr1)].replace('.', '')) - 1
    rec.append(str(pos))     # CDR1 end
    pos += 1
    rec.append(str(pos))     # FWR2 start
    pos += len(seq[slice(*imgt_fr2)].replace('.', '')) - 1
    rec.append(str(pos))     # FWR2 end
    pos += 1
    rec.append(str(pos))     # CDR2 start
    pos += len(seq[slice(*imgt_cdr2)].replace('.', '')) - 1
    rec.append(str(pos))     # CDR2 end
    pos += 1
    rec.append(str(pos))     # FWR3 start
    pos += len(seq[slice(*imgt_fr3)].replace('.', '')) - 1
    rec.append(str(pos))     # FWR3 end
    rec.append(chain)    # chain type
    rec.append('0')       # coding frame start
    return '\t'.join(rec) + '\n'
    


def from_json(args):
    with open(args.ref_file, 'r') as fo:
        ref = json.load(fo)

    if isinstance(ref["GermlineSet"], list):
        germline_set = ref["GermlineSet"][0]
    else:
        germline_set = ref["GermlineSet"]

    with open(args.ndm_file, 'w') as fo:
        for allele in germline_set['allele_descriptions']:
            if allele['sequence_type'] == 'V':
                if 'v_gene_delineations' not in allele:
                    print(f"Omitting allele {allele['label']} as no delineations are specified")
                    continue

                for delineation in allele['v_gene_delineations']:
                    if delineation['delineation_scheme'] == 'IMGT':
                        for coord in ['fwr1_start', 'fwr1_end', 'cdr1_start', 'fwr2_end', 'cdr2_start', 'fwr3_end']:
                            if coord not in delineation or delineation[coord] is None:
                                print(f"Omitting incomplete allele {allele['label']}")
                                continue

                        rec = record_from_json(allele['label'], delineation, args.chain)
                        
                        # compare with gapped sequence derivation and warn if there is a mismatch
                        # can't make this comparison any longer because we don't know the coordinates used in the gapped alignment
                        # rec_from_gapped = record_from_gapped_seq(allele['label'], delineation['aligned_sequence'], args.chain)

                        # if rec != rec_from_gapped:
                        #     print(f"Warning: mismatch between IMGT gapped sequence and CDR delineation for allele {allele['label']}")

                        fo.write(rec)

    


def record_from_json(label, d, chain):
    return(f"{label}\t{d['fwr1_start']}\t{d['fwr1_end']}\t{d['cdr1_start']}\t{d['cdr1_end']}\t{d['fwr2_start']}\t{d['fwr2_end']}\t{d['cdr2_start']}\t{d['cdr2_end']}\t{d['fwr3_start']}\t{d['fwr3_end']}\t{chain}\t0\n")


if __name__ == '__main__':
    main()  
