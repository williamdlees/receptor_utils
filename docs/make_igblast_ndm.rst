.. _make_igblast_ndm_label:

make_igblast_ndm
================

.. argparse::
   :filename: ../src/scripts/make_igblast_ndm.py
   :func: get_parser
   :prog: make_igblast_ndm

Notes
=====

If the reference set is available in AIRR-C JSON format, it is recommended to use this format as it will ensure correct placement of the CDRs using coordinates for each sequence from the file. 
All germline sets on `OGRDB <https://ogrdb.airr-community.org>`_ are available in this format.

If the fasta format is used, V sequences must be IMGT gapped. By default, make_igblast_ndm assumes that the gapped alignment follows the  `IMGT numbering scheme <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumbering.html>`_ 
with no inserted codons. Please note that some reference sets available from IMGT, such as those for rhesus macaque, do contain insertions in the alignment. In this case, you can use the ``--cdr_coords`` option to 
specify coordinates of the CDRs. If in doubt, please analyse one sequence from the gapped set in IMGT V-QUEST and review the codon alignment: if there are insertions, denoted by a codon number followed by A, B, C, etc.,
you will need to specify the CDR coordinates using the ``--cdr_coords`` option. An alternative is to review the sequences in the reference set and remove the insertions. This may be an option, depending upon
the application, the number of germline sequences that have the insertion, and their usage in your samples. We offer a specific tool, :ref:`fix_macaque_gaps_label` to conduct this realignment for IMGT rhesus macaque sets.

Please note that downstream analysis tools may assume that IMGT-gapped sequences follow the default alignment with no inserted codons. If you use a tool that requires a gapped set of V sequences and there is no 
provision to specify CDR coordinates, this is likely to be the case.
