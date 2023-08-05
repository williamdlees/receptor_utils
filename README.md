# receptor_utils

Some tools I find useful for working with Ig receptor sequences, including support for allele sequence naming and 
the creation of custom IgBlast databases. Please see the [documentation](https://williamdlees.github.io/receptor_utils/_build/html/introduction.html)
for further details.

Changes in version 0.0.38:
annotate_j and make_igblast_ndm will now accept a germline set in AIRR Community JSON format, as an alternative to providing the set in FASTA format.

Changes in version 0.0.39:
Improve reporting of issues with conserved residues
Change URL for fetching IMGT reference sets to use https

Changes in version 0.0.40:
The submodule name receptor_utils.number_ighv has been changed to receptor_utils.number_v to reflect its wider scope. The old name will continue
to work for the time being but will raise a deprecation warning.

in receptor_utils.simple_bio_seq, write_fasta(seqs, filename) has become write_fasta(filename, seqs) for consistency with write_csv. The old
calling pattern will continue to work for the time being but will raise a deprecation warning.



