.. _download_germline_set:

download_germline_set
=====================

This utility will download reference sequences from the Open Germline Receptor Database (`OGRDB <https://ogrdb.airr-community.org>`_).

.. argparse::
   :filename: ../src/scripts/download_germline_set.py
   :func: get_parser
   :prog: download_germline_set

Format Options
--------------

* JSON format: AIRR-C JSON format (default)
* Single FASTA file: all V(D)J sequences in a single FASTA file
* Multiple FASTA files: V, D, J, and gapped V sequences in separate FASTA files
* IgBLAST format: Multiple fasta files, plus IgBLAST germline configuration files

By default, germline sets are downloaded in `AIRR-C JSON format <https://docs.airr-community.org/en/stable/datarep/germline.html/>`_. The advantage of this format
is that it contains full information on the germline set including delineation of the V-sequence CDRs and delieation of the J-sequences. This means that it can,
in principle, be loaded into an annotation tool without needing any additional information to be provided to the tool. The json format is also recommended for use 
with :ref:`annotate_j_label` and :ref:`make_igblast_ndm_label` utilities, which create germline configuration files for IgBLAST. The AIRR-C JSON format can also be 
specified with the argument -f AIRRC-JSON.

Another option is to download all V(D)J sequences into a single FASTA file. This can be specified by the -f SINGLE-FG argument, which will download gapped V sequences, or
-f SINGLE-FU, which will download ungapped V sequences. The filename, by default, is, respectively, <species>_<locus>_gapped.fasta, or <species>_<locus>.fasta. 
The -p argument can be used to replace the <species>_<locus>_ prefix with a custom prefix.

The sequences can also be downloaded into four FASTA files. These are named, by default, <species>_<locus>_V.fasta, <species>_<locus>_D.fasta, <species>_<locus>_J.fasta, 
<species>_<locus>_V_gapped.fasta. This is specified by the -f MULTI_F argument. The -p argument can be used to replace the <species>_<locus>_ prefix with a custom prefix.

Finally, the -f MULTI-IGBLAST option will download the four multi files as above, and also create .ndm and .aux files for use with IgBLAST. For further details 
on use with IgBLAST please see :ref:`airrc_sets_with_igblast`.

Examples
--------

Download the latest version of the human IGK germline set in AIRR-C JSON format::

   download_germline_set -"Homo sapiens" IGK

Download the latest version of the human IGH germline set in single FASTA format, with gapped V sequences::

   download_germline_set "Homo sapiens" IGH -f SINGLE-FG

Download the latest version of the mouse C57BL/6 IGH germline set in multiple FASTA format::

   download_germline_set "Mus musculus" IGH -n "C57BL/6 IGH" -f MULTI-F

Note that the utility will try to provide helpful information in the event of a command error::

   >download_germline_set "Mus musculus" IGH -n "C57BL/6" -f MULTI-F
   https://ogrdb.airr-community.org/api_v2/germline/species
   Mus musculus: 10090
   Error: set C57BL/6 not found for species 10090 locus IGH. Available sets:
   BALB/c IGHV, C57BL/6 IGH, C57BL/6 IGHV, CAST/EiJ IGH, LEWES/EiJ IGH, MSM/MsJ IGH, NOD/ShiLtJ IGH, PWD/PhJ IGH, BALB/c IGH

Download human IGH files for IgBLAST, using the custom filename prefix AIRRC_IGH::

   download_germline_set "Homo sapiens" IGH -f MULTI-IGBLAST -p AIRRC_IGH