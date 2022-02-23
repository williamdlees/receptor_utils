# Provide a name for a novel allele sequence, using the VDJbase approach as outlined
# at https://wordpress.vdjbase.org/index.php/vdjbase_help/airr-seq-data-allele-names/
# For V sequences, the sequence must be full-length at the 5' end, or IMGT-gapped.

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE


from receptor_utils import simple_bio_seq as simple
from receptor_utils import number_ighv

# Make a name for the novel allele, given its gapped or ungapped sequence
# returns a tuple (name, gapped_sequence)
#
# novel_seq: the full-length sequence to name, which may be either gapped or ungapped in the case of a V-gene
# ref_set: dict of reference genes (gapped in the case of v-genes), in the format returned by read_fasta
# v_gene: True if we are naming a v_gene
# gaps in V-sequences are '.', as in the IMGT reference set


def name_novel(novel_seq, ref_set, v_gene=True):
    novel_seq = novel_seq.upper()
    notes = ''

    # gap the sequence (even if already gapped) so that we get the notes

    if v_gene:
        if '.' in novel_seq:
            novel_seq = novel_seq.replace('.', '')

        ungapped_ref_set = {}
        for k, v in ref_set.items():
            ungapped_ref_set[k] = v.replace('.', '')
        (novel_seq, aa, notes) = number_ighv.gap_sequence(novel_seq, ref_set, ungapped_ref_set)

    # Check for truncation at 5' end and construct ambiguous reference set if necessary

    if novel_seq[0] == '.':
        for i in range(0, len(novel_seq)):
            if novel_seq[i] != '.':
                break
        ref_set = build_ambiguous_ref(ref_set, i)

    closest_ref_name = simple.closest_ref(novel_seq, ref_set)[0]
    closest_ref_seq = ref_set[closest_ref_name]

    if len(closest_ref_seq) < len(novel_seq):
        closest_ref_seq += '.' * (len(novel_seq) - len(closest_ref_seq))

    diffs = {}
    class Diff:
        pass

    for i in range(0, len(novel_seq)):
        if closest_ref_seq[i] != novel_seq[i]:
            d = Diff()
            d.pos = i
            d.ref = closest_ref_seq[i]
            d.novel = novel_seq[i]
            diffs[i] = d

    # trim any diffs that are past the end of the novel sequence

    for i in range(len(novel_seq) - 1, 0):
        if novel_seq[i] != '.':
            break
    novel_end = i

    for i in list(diffs.keys()):
        if i > novel_end:
            del diffs[i]

    # add diffs if the closest ref is longer than the novel seq

    if len(closest_ref_seq) > novel_end + 1:
        for j in range(novel_end, len(closest_ref_seq) - 1):
            d = Diff()
            d.pos = j
            d.ref = closest_ref_seq[j]
            d.novel = '-'
            diffs[j] = d


    # consolidate any adjacent diffs

    prev_i = None
    for i, d in list(diffs.items()):
        if prev_i and diffs[prev_i].pos + len(diffs[prev_i].novel) == i:
            diffs[prev_i].ref += closest_ref_seq[i]
            diffs[prev_i].novel += novel_seq[i] if i < len(novel_seq) - 1 else '-'
            del diffs[i]
        else:
            prev_i = i

    # split diffs that are too short to consolidate

    for i, d in list(diffs.items()):
        if 1 < len(d.novel) < 4 and '.' not in d.novel and '.' not in d.ref:
            for j in range(1, len(d.novel)):
                diffs[i + j] = Diff()
                diffs[i + j].pos = d.pos + j
                diffs[i + j].ref = d.ref[j]
                diffs[i + j].novel = d.novel[j]
            d.ref = d.ref[0]
            d.novel = d.novel[0]

    # form the suffix

    suffs = []
    for i in sorted(diffs.keys()):
        d = diffs[i]
        if len(d.novel) == 1 and d.ref != '.' and d.novel != '.':
            suffs.append('%s%d%s' % (d.ref.lower(), d.pos + 1, d.novel.lower()))
        else:
            suffs.append('%d%s%d' % (d.pos + 1, d.novel.lower(), d.pos + len(d.novel)))

    novel_name = closest_ref_name + '_' + '_'.join(suffs) if len(suffs) else closest_ref_name
    return(novel_name, novel_seq, notes)


# Build (and name) a reference set for sequences truncated at the 5' end
# Start is 0-based index of first available bp


def build_ambiguous_ref(full_ref, start):
    rev_ref = {}
    for name, seq in full_ref.items():
        seq = '.' * start + seq[start:]
        if seq not in rev_ref:
            rev_ref[seq] = []
        rev_ref[seq].append(name)

    ref = {}
    for seq, names in rev_ref.items():
        if'IGHV3-15*08' in names:
            print('foo')

        if len(names) == 1:
            ref[names[0]] = seq
        else:
            genes = {}
            for name in names:
                gene, allele = name.split('*')
                if gene not in genes:
                    genes[gene] = []
                genes[gene].append(allele)
            if len(genes) == 1:
                name = gene + '*' + '_'.join(sorted(genes[gene]))
                ref[name] = seq
            else:
                gene_list = sorted(genes.keys())
                name = gene_list[0] + '*' + '_'.join(sorted(genes[gene_list[0]]))
                fam_0 = gene_list[0].split('-')[0]
                for i in range(1, len(gene_list)):
                    name += '_'
                    fam = gene_list[i].split('-')[0]
                    num = gene_list[i].replace(fam + '-', '')
                    if fam == fam_0:
                        alleles = [num + '.' + a for a in genes[gene_list[i]]]
                        name += '_'.join(alleles)
                    else:
                        fam = gene_list[i][4:]
                        alleles = [fam + '.' + a for a in genes[gene_list[i]]]
                        name += '_'.join(alleles)
                ref[name] = seq

    return ref


def run_tests():
    test_novels = simple.read_fasta('tests/IGHV3-20.04.fasta')
    ref_set = simple.read_fasta('refs/Homo_Sapiens_IGHV_gapped.fasta')

    for name, seq in test_novels.items():
        novel_name = name_novel(seq, ref_set)
        if novel_name != name:
            print('Error: %s != %s' % (novel_name, name))


if __name__ == "__main__":
    run_tests()

