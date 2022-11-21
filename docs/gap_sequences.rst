.. _gap_sequences_label:

gap_sequences
=============

Gap the sequences in a FASTA file, using the closest sequence in a gapped reference set as template.

The tool will only use as a template those sequences in the reference set that are complete at the 5' end and extend at least to the
second cysteine. Any others will not be used as a template, and neither will any that are apparently non-functional.
Warnings will be given of any sequences with these issues. In most circumstances (that is, if they are only a small
fraction of the total sequences in the reference set), these warnings can be ignored.

.. argparse::
   :filename: ../src/scripts/gap_inferred.py
   :func: get_parser
   :prog: gap_sequences
