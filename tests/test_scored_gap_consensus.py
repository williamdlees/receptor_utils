from receptor_utils import simple_bio_seq as simple

seqs = [
    'abcdef',
    'abcded',
    'aacdef',
]

seq_dict = {}
for i in range(len(seqs)):
    seq_dict[str(i)] = seqs[i]

c, s = simple.scored_consensus(seq_dict)

breakpoint()