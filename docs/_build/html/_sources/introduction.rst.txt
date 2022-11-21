Introduction
============

This package contains utilities and code which I have found useful in working with IG/TR receptor sequences.

Installation
------------

``pip install receptor-utils``

The module requires `Biopython <https://biopython.org/>`_, which needs to be installed separately.


Naming and aligning V, D and J sequences
----------------------------------------

Particularly when dealing with novel inferences of germline genes, but in some other contexts as well, when given an expressed or genomic sequence, it is useful to know which germline V, D or J sequence it
is closest to, and what SNPs it contains. An allied task is to align such a sequence using the `IMGT alignment <https://www.imgt.org/IMGTindex/AlleleAlignments.php>`_. 

:ref:`name_allele_label` provides a name for a given sequence. The name includes the nearest germline allele (taken from a reference set supplied to the tool) and incorporating any SNPs. The naming scheme that it
uses is described `here <https://wordpress.vdjbase.org/index.php/vdjbase_help/airr-seq-data-allele-names/>`_.

:ref:`gap_sequences_label` will gap a supplied set of sequences according to the IMGT alignment, warning of any missing conserved residues. It employs a gapped reference set, which must be
provided. For each sequence to be gapped, the tool identifies the closest reference sequence, and uses it as a template.


Using custom databases with IgBlast
-----------------------------------

The package contains some tools to support this: they are covered in :ref:`the next section <custom_igblast_label>`.

Convenience Tools
-----------------

There are some smaller 'convenience' tools that may come in useful, for example to :ref:`extract sequences from an IMGT reference set <extract_refs_label>`
and to :ref:`reverse-complement an nt sequence <rev_comp_label>`.
They are nothing special but may save some time. All tools are listed under :ref:`command_line_label`.

simple_bio_seq API
------------------

This is a set of simple wrappers around commonly-used functions in Biopython which suit my use case and are used in the tools above. They support, for example,
one-line reading and writing of FASTA files and comma/tab-separated files, and manage sequences as strings in a dict, which I find simplifies code compared the
Biopython functions, at the expense of flexibility which I rarely need. The are documented under :ref:`API <api_label>`.