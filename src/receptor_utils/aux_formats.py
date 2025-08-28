# Functions to create proprietary files for aligners (IgBLAST, MiXCR) from reference data

import json
import re
from collections import namedtuple
from receptor_utils import simple_bio_seq as simple


def ndm_from_json(ref, ndm_file, chain, delineation_scheme='IMGT'):
    if isinstance(ref["GermlineSet"], list):
        germline_set = ref["GermlineSet"][0]
    else:
        germline_set = ref["GermlineSet"]

    if 'items' in germline_set:
        germline_set = germline_set['items'][0]

    ret = False
    with open(ndm_file, 'w') as fo:
        for allele in germline_set['allele_descriptions']:
            if allele['sequence_type'] == 'V':
                if 'v_gene_delineations' not in allele:
                    print(f"Omitting allele {allele['label']} as no delineations are specified")
                    continue

                if not ret:
                    fo.write('#name\tFWR1 start\tFWR1 stop\tCDR1 start\tCDR1 stop\tFWR2 start\tFWR2 stop\tCDR2 start\tCDR2 stop\tFWR3 start\tFWR3 stop\tchain type\tframe start\n')
                    ret = True

                for delineation in allele['v_gene_delineations']:
                    if delineation['delineation_scheme'] == delineation_scheme:
                        for coord in ['fwr1_start', 'fwr1_end', 'cdr1_start', 'fwr2_end', 'cdr2_start', 'fwr3_end']:
                            if coord not in delineation or delineation[coord] is None:
                                print(f"Omitting incomplete allele {allele['label']}")
                                continue

                        rec = record_from_json(allele['label'], delineation, chain)

                        # compare with gapped sequence derivation and warn if there is a mismatch
                        # can't make this comparison any longer because we don't know the coordinates used in the gapped alignment
                        # rec_from_gapped = record_from_gapped_seq(allele['label'], delineation['aligned_sequence'], args.chain)

                        # if rec != rec_from_gapped:
                        #     print(f"Warning: mismatch between IMGT gapped sequence and CDR delineation for allele {allele['label']}")

                        fo.write(rec)

    return ret


def record_from_json(label, d, chain):
    return f"{label}\t{d['fwr1_start']}\t{d['fwr1_end']}\t{d['cdr1_start']}\t{d['cdr1_end']}\t{d['fwr2_start']}\t{d['fwr2_end']}\t{d['cdr2_start']}\t{d['cdr2_end']}\t{d['fwr3_start']}\t{d['fwr3_end']}\t{chain}\t0\n"


def ndm_from_fasta(ref_file, ndm_file, chain, cdr_coords):
    seqs = simple.read_fasta(ref_file)
    # desired headers:
    # gene/allele name, FWR1 start, FWR1 stop, CDR1 start, CDR1 stop, FWR2 start, FWR2 stop, CDR2 start, CDR2 stop, FWR3 start, FWR3 stop, chain type, coding frame start.
    # FWR/CDR positions are 1-based while the coding frame start positions are 0-based

    omissions = False

    with open(ndm_file, 'w') as fo:
        fo.write('#name\tFWR1 start\tFWR1 stop\tCDR1 start\tCDR1 stop\tFWR2 start\tFWR2 stop\tCDR2 start\tCDR2 stop\tFWR3 start\tFWR3 stop\tchain type\tframe start\n')

        for label, seq in seqs.items():

            # If the sequence is truncated at the 5' end, don't include the sequence

            if seq[0] == '.':
                print(f"Omitting allele {label} as it is truncated at the 5' end")
                omissions = True
                continue

            fo.write(record_from_gapped_seq(label, seq, chain, cdr_coords))

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


def aux_from_json(ref, out_file, verbose):
    if isinstance(ref["GermlineSet"], list):
        germline_set = ref["GermlineSet"][0]
    else:
        germline_set = ref["GermlineSet"]

    if 'items' in germline_set:
        germline_set = germline_set['items'][0]

    annotations = []

    for allele in germline_set['allele_descriptions']:
        if allele['sequence_type'] == 'J':

            try:
                gene_start = int(allele['gene_start'])
                gene_end = int(allele['gene_end'])
                j_cdr3_end = int(allele['j_cdr3_end'])
                j_codon_frame = int(allele['j_codon_frame'])    
            except Exception as e:
                print(f"Unrecognisable parameter in JSON file for allele {allele['label']}: {e}. Skipping.")
                continue

            if (j_cdr3_end - j_codon_frame) % 3 != 0:
                print(f"CDR3 end for {allele['label']} is not on a codon boundary. Skipping.")
                continue

            annotations.append({
                '#name': allele['label'],
                'j_codon_frame': allele['j_codon_frame'] - 1,
                'chain_type': f"{allele['label'][3]}{allele['label'][2]}",
                'j_cdr3_end': j_cdr3_end - 2,
                'extra_bps': (gene_end - gene_start - j_codon_frame - 1 ) % 3,
            })

    annotations = sorted(annotations, key=lambda x: x['#name'])

    ret = False
    if annotations:
        simple.write_csv(out_file, annotations, delimiter='\t')
        ret = True

    return ret


def aux_from_fasta(seq_file, out_file, verbose):
    seqs = simple.read_fasta(seq_file)
    aux_from_seqs(seqs, out_file, verbose)


def aux_from_seqs(seqs, out_file, verbose):
    if verbose:
        print('Annotations use MiAIRR definitions for codon_frame, j_cdr3_end')

    annotations = []

    for seq_name, seq in seqs.items():
        solutions = []
        Result = namedtuple('Result', ['nt_seq', 'frame', 'cdr3_pos', 'trans', 'extra_bps', 'stop_codon', 'canonical_junc', 'motif'])
        for frame in range(3):
            trans = simple.translate(seq[frame:])

            # if '*' in trans:      allow stop codon as may be excised during junction formation
            #    continue

            for m in re.finditer('[WF]G.G', trans):
                solutions.append(Result._make([
                    seq,
                    frame+1,
                    frame + 3*m.start() + 1,
                    trans,
                    (len(seq) - frame) % 3,
                    '*' in trans,
                    True,
                    m.group(0)
                ]))

            if not solutions:
                for m in re.finditer('[CVWF]G.G', trans):     # found in TRJP1, TRJP2
                    solutions.append(Result._make([
                        seq,
                        frame+1,
                        frame + 3*m.start() + 1,
                        trans,
                        (len(seq) - frame) % 3,
                        '*' in trans,
                        False,
                        m.group(0)
                    ]))
                    #print(f'{seq_name}: solution with non-standard motif [CVWF][AG].G')

            if not solutions:
                for m in re.finditer('[WF][AG].G', trans):     # found in TRJP1, TRJP2
                    solutions.append(Result._make([
                        seq,
                        frame+1,
                        frame + 3*m.start() + 1,
                        trans,
                        (len(seq) - frame) % 3,
                        '*' in trans,
                        False,
                        m.group(0)
                    ]))
                    #print(f'{seq_name}: solution with non-standard motif [CWF]G.G')

            if not solutions:
                for m in re.finditer('[WF]G.[NG]', trans):     # found in TRJP1, TRJP2
                    solutions.append(Result._make([
                        seq,
                        frame+1,
                        frame + 3*m.start() + 1,
                        trans,
                        (len(seq) - frame) % 3,
                        '*' in trans,
                        False,
                        m.group(0)
                    ]))
                    #print(f'{seq_name}: solution with non-standard motif [WF]G.[NG]')

        if solutions:
            if len(solutions) > 1:
                # if there is a solution with a canonical junction, prefer it
                preferred = [s for s in solutions if s.canonical_junc]
                if preferred:
                    solutions = preferred
                if len(solutions) > 1:
                    # if there is a solution without a stop codon, take it
                    preferred = [s for s in solutions if not s.stop_codon]
                    if preferred:
                        solutions = preferred
                    if len(solutions) > 1:
                        print(f'{seq_name}: multiple solutions found, selecting one')
                        solutions = [solutions[0]]
                
            for solution in solutions:
                if '*' in solution.trans:
                    print(f'{seq_name}: solution with stop codon in translation')

                if not solution.canonical_junc:
                    print(f'{seq_name}: solution with non-canonical motif {solution.motif}')

                if verbose:
                    print(f'{seq_name}: j_codon_frame {solution.frame} j_cdr3_end {solution.cdr3_pos} () exta bps {solution.extra_bps}')

                    loc_count = '      '
                    for num in range(len(solution.trans)):
                        loc_count += f' {num+1:<2}'
                    print(loc_count)

                    aa = '      '
                    for num in range(len(solution.trans)):
                        aa += f' {solution.trans[num]} '
                    print(aa)

                    nt = f'  {solution.nt_seq[:solution.frame-1]:>3} {solution.nt_seq[solution.frame-1:]}'
                    print(nt + '\n\n')

                annotations.append({
                    '#name': seq_name,
                    'j_codon_frame': solution.frame - 1,
                    'chain_type': f'{seq_name[3]}{seq_name[2]}',
                    'j_cdr3_end': solution.cdr3_pos-2,
                    'extra_bps': solution.extra_bps,
                })
        else:
            print(f'{seq_name}: no solutions')

    simple.write_csv(out_file, annotations, delimiter='\t')
