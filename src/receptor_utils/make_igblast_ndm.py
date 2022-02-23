#!python
# Make igblast ndm file from a gapped reference set

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


import csv
from receptor_utils import simple_bio_seq as simple
import argparse


def main():
    parser = argparse.ArgumentParser(description='Make igblast ndm file from an IMGT-gapped variable-region reference set')
    parser.add_argument('ref_file', help='gapped reference set (fasta)')
    parser.add_argument('chain', help='chain as required by igblast (eg VH)')
    parser.add_argument('ndm_file', help='output file (ndm)')
    args = parser.parse_args()
    
    seqs = simple.read_fasta(args.ref_file)
    
    # Python-based ranges of IMGT elements
    # <---------------------------------- FR1-IMGT -------------------------------->______________ CDR1-IMGT ___________<-------------------- FR2-IMGT ------------------->___________ CDR2-IMGT ________<----------------------------------------------------- FR3-IMGT ----------------------------------------------------> CDR3-IMGT
    imgt_fr1 = (0, 78)
    imgt_cdr1 = (78, 114)
    imgt_fr2 = (114, 165)
    imgt_cdr2 = (165, 195)
    imgt_fr3 = (195, 312)
    
    # desired headers:
    # gene/allele name, FWR1 start, FWR1 stop, CDR1 start, CDR1 stop, FWR2 start, FWR2 stop, CDR2 start, CDR2 stop, FWR3 start, FWR3 stop, chain type, coding frame start.
    # FWR/CDR positions are 1-based while the coding frame start positions are 0-based
    
    with open(args.ndm_file, 'w') as fo:
        for id, seq in seqs.items():
            rec = []
            rec.append(id)
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
            rec.append(args.chain)    # chain type
            rec.append('0')       # coding frame start
            fo.write('\t'.join(rec) + '\n')
