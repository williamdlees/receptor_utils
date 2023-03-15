.. _merge_fasta_label:

merge_fasta
===========

This utility will merge sequences from two fasta files. The output will contain all sequences from file1, and any sequences
from file2 that are not present in file1 (the sequences themselves are compared, rather than the sequence ids).

.. argparse::
   :filename: ../src/scripts/merge_fasta.py
   :func: get_parser
   :prog: merge_fasta
