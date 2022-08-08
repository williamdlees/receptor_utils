# Wrapper around BioSeq functions for simple applications
# principles:
# - store sequences as strings, use dicts for collectiopns
# - all sequences are coerced to upper case on input
# - iterators are coreced into lists for ease of debugging
# basically just make things simple for cases where we don't need to do more

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
import random
import csv
from collections import Counter


def chunks(l, n):
    """ Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


def read_fasta(infile):
    res = {}
    recs = SeqIO.parse(infile, 'fasta')
    for rec in recs:
        res[rec.id] = str(rec.seq).upper()
    return res


def write_fasta(seqs, outfile):
    SeqIO.write(toSeqRecords(seqs), outfile, 'fasta')


def read_single_fasta(infile):
    res = ''
    recs = SeqIO.parse(infile, 'fasta')
    for rec in recs:
        return str(rec.seq).upper()


# read V,D,J regions from one or more species from an IMGT reference file
# species should be a ***list*** of species
# returns a nested dict: [species][chain]
def read_imgt_fasta(infile, species, chains=('IGHV', 'IGHD', 'IGHJ', 'CH'), functional_only=False):
    res = {}
    for sp in species:
        res[sp] = {}
        for chain in chains:
            res[sp][chain] = {}

    recs = SeqIO.parse(infile, 'fasta')

    for rec in recs:
        for sp in species:
            for chain in chains:
                if ('-REGION' in rec.description or chain == 'CH' or 'EX' in rec.description) and sp in rec.description and chain in rec.description:
                    if not functional_only or '|F|' in rec.description:
                        if 'V' not in chain or len(rec.seq) > 280:       # strip out obviously incomplete Vs
                            name = rec.description.split('|')[1]
                            type = rec.description.split('|')[4]

                            if chain == 'CH' or ('EX' in rec.description and chain in rec.description):
                                gene, allele = name.split('*')
                                res[sp][chain][gene + '_' + type + '*' + allele] = str(rec.seq).upper()
                            else:
                                res[sp][chain][name] = str(rec.seq).upper()
    return res


def sample_fasta(seqs, number):
    keys = seqs.keys()

    if len(keys) <= number:
        return seqs

    sample = random.sample(keys, number)
    sample_seqs = {}
    for key in sample:
        sample_seqs[key] = seqs[key]
    return sample_seqs


def translate(seq, truncate=True, ignore_partial_codon=True):
    seq_len = len(seq)
    residue = seq_len % 3

    if residue != 0:
        if truncate:
            seq = seq[:seq_len - residue]
        else:
            seq += 'N' * (3 - residue)
            
    if ignore_partial_codon:
        codons = [codon for codon in chunks(seq, 3)]
        i = 0

        for i in range(len(codons)):
            if '-' in codons[i] or '.' in codons[i]:
                codons[i] = '---'
            
        seq = ''.join(codons)

    return str(Seq(seq).translate())


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def nt_diff(s1, s2):
    diffs = 0
    row = list(zip(s1.upper(), s2.upper()))

    for p in range(len(row)):
        if row[p][0] != row[p][1] and (row[p][0] not in ['X', 'N'] and row[p][1] not in ['X', 'N']):
            diffs += 1

    diffs += abs(len(s1) - len(s2))
    return diffs


# convert our dict format to SeqRecords
def toSeqRecords(seqs):
    recs = []
    for name, seq in seqs.items():
        recs.append(SeqRecord(Seq(seq), id=name, description=''))
    return recs


def dumb_consensus(seqs, threshold=0.7):
    msa = MultipleSeqAlignment(toSeqRecords(seqs))
    summary_align = AlignInfo.SummaryInfo(msa)
    return str(summary_align.gap_consensus(threshold=threshold))


def scored_consensus(seqs, threshold=0.7):
    msa = MultipleSeqAlignment(toSeqRecords(seqs))
    summary_align = AlignInfo.SummaryInfo(msa)
    return scored_gap_consensus(msa, threshold=threshold)


# Adapted from BioPython.Bio.ALign
# return min_score: this is the lowest proportion of reads that 'won' at any position
def scored_gap_consensus(alignment, threshold=0.7, ambiguous="X", require_multiple=False):
    """Output a fast consensus sequence of the alignment, allowing gaps.
    Same as dumb_consensus(), but allows gap on the output.
    Things to do:
     - Let the user define that with only one gap, the result
       character in consensus is gap.
     - Let the user select gap character, now
       it takes the same as input.
    """
    consensus = ""

    # find the length of the consensus we are creating
    con_len = alignment.get_alignment_length()
    con_rows = len(alignment)

    # go through each seq item
    min_score = 1.0
    for n in range(con_len):
        # keep track of the counts of the different atoms we get
        atom_dict = Counter()
        num_atoms = 0

        for record in alignment:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            try:
                c = record[n]
            except IndexError:
                continue
            atom_dict[c] += 1

            num_atoms += 1

        max_atoms = []
        max_size = 0

        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)

        if require_multiple and num_atoms == 1:
            consensus += ambiguous
            min_score = min(min_score, 1/con_rows)
        elif (len(max_atoms) == 1) and (
            (float(max_size) / float(num_atoms)) >= threshold
        ):
            consensus += max_atoms[0]
            min_score = min(min_score, float(max_size) / float(con_rows))
        else:
            consensus += ambiguous
            min_score = min(min_score, float(max_size) / float(con_rows))

    return (Seq(consensus), min_score)


# Find the closest reference sequence to a query sequence
# Returns a list of equally close reference names
def closest_ref(seq, ref):
    closest_diff = 999
    closest_names = []
    for ref_name, ref_seq in ref.items():
        diff = nt_diff(seq, ref_seq)
        if diff < closest_diff:
            closest_diff = diff
            closest_names = [ref_name]
        elif diff == closest_diff:
            closest_names.append(ref_name)
    return closest_names


# Read csv file into a list of dicts
def read_csv(file, delimiter=None):
    ret = []
    with open(file, 'r') as fi:
        if delimiter:
            reader = csv.DictReader(fi, delimiter=delimiter)
        else:
            reader = csv.DictReader(fi)
        for row in reader:
            ret.append(row)

    return ret


# Write csv file given a list of dicts. Fieldnames are taken from the first row
def write_csv(file, rows, delimiter=None):
    if not rows:
        return

    fieldnames = rows[0].keys()
    with open(file, 'w', newline='') as fo:
        if delimiter:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter=delimiter)
        else:
            writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()
        
        for row in rows:
            writer.writerow(row)
            