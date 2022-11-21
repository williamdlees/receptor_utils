# Copyright (c) 2022 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE

# Name an allele by comparing to a reference set. Construct a novel name if necessary

from receptor_utils import simple_bio_seq as simple
from receptor_utils import novel_allele_name
import argparse


def find_gene(name):
    gene = name
    if '*' in name:
        gene = gene.split('*')[0]
    return gene


def get_parser():
    parser = argparse.ArgumentParser(description='Name an allele by comparing to a reference set. Construct a novel name containing SNPs if necessary')
    parser.add_argument('ref_set', help='reference set (FASTA) (V genes should be IMGT-gapped)')
    parser.add_argument('sequence', help='input sequence')
    parser.add_argument('-g', '--gene', metavar='GENE_NAME', help='only consider reference sequences from the specified gene', dest='gene')
    parser.add_argument('-r', '--rev_comp', help='reverse-complement sequence before processing', action='store_true', dest='rev_comp')
    return parser


def main():
    args = get_parser().parse_args()

    gene_refs = simple.read_fasta(args.ref_set)
    
    if args.gene:
        gene_refs = {k: v for k, v in gene_refs.items() if find_gene(k) == args.gene}
        
    if len(gene_refs) == 0:
        print(f"No reference genes to compare")
        quit()
        
    seq = args.sequence
    if args.rev_comp:
        seq = simple.reverse_complement(seq)

    novel_name, novel_seq, notes = novel_allele_name.name_novel(seq, gene_refs, len(seq) > 60)
           
    print(novel_name)
    
    if '.' in novel_seq:
        print(f"Gapped sequence:\n{novel_seq}")

    if notes:
        print(f"Notes:\n{notes}")
