.. _download_germline_set:

download_germline_set
=====================

This utility will download reference sequences from the Open Germline Receptor Database (OGRDB) at https://ogrdb.airr-community.org.

.. argparse::
   :filename: ../src/scripts/download_germline_set.py
   :func: get_parser
   :prog: download_germline_set

Explanation of the -f argument
==============================

The utility has been designed to download germline sequences in a flexible manner.

By default, germline sets are downloaded in `AIRR-C JSON format <https://docs.airr-community.org/en/stable/datarep/germline.html/>`_. The advantage of this format
is that it contains full information on the germline set including delineation of the V-sequence CDRs and delieation of the J-sequences. This means that it can,
in principle, be loaded into an annotation tool without needing any additional information to be provided to the tool. The json format is also recommended for use 
with :ref:`annotate_j_label` and :ref:`make_igblast_ndm_label` utilities, which create germline configuration files for IgBLAST. The AIRR-C JSON format can also be 
specified with the argument -f AIRRC-JSON.

Another option is to download all V(D)J sequences into a single FASTA file. This can be specified by the -f SINGLE-FG argument, which will download gapped V sequences, or
-f SINGLE-FU, which will download ungapped V sequences. The filename, by default, is, respectively, <species>_<locus>_gapped.fasta, or <species>_<locus>.fasta. 
The -p argument can be used to replace the <species>_<locus>_ prefix with a custom prefix.

The final option is to download the sequences into four FASTA files. These are named, by default, <species>_<locus>_V.fasta, <species>_<locus>_D.fasta, <species>_<locus>_J.fasta, 
<species>_<locus>_V_gapped.fasta. This is specified by the -f MULTI_F argument. The -p argument can be used to replace the <species>_<locus>_ prefix with a custom prefix.
