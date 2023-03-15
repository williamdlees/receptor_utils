.. _fix_macaque_gaps_label:

fix_macaque_gaps
================

IMGT's gapped FASTA sequences for rhesus macaque IGHV, IGLV and IGKV contain inserted codons relative to the sequences for
most other species (I am not aware of any other species for which this is the case, but have not conducted an
exhaustive search). These additional codons in macaque sequences cause problems for downstream tools such as Change-o
that expect features such as CDRs and conserved residues to be present at specific co-ordinates. This utility restores
the 'standard' alignment by removing the inserted codons. A small number of sequences in the IMGT sets use these
codons: the tool lists those sequences and removes them from the set. Personally I have not observed these sequences
in rhesus repertoires, but this is something that you may wish to check before using the set produced by the tool.

.. argparse::
   :filename: ../src/scripts/fix_macaque_gaps.py
   :func: get_parser
   :prog: fix_macaque_gaps
