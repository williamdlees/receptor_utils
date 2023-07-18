Introduction
============

This package contains utilities and code which are useful in working with IG/TR receptor sequences.

Installation
------------

``pip install receptor-utils``

The module requires `Biopython <https://biopython.org/>`_, which needs to be installed separately.

Naming and aligning V, D and J sequences
----------------------------------------

:ref:`name_allele_label` - Provide a name for a given sequence. The name includes the nearest
germline allele (taken from a reference set supplied to the tool) and
incorporating any SNPs. The naming scheme that it uses is described
`here <https://wordpress.vdjbase.org/index.php/vdjbase_help/airr-seq-data-allele-names/>`_

:ref:`gap_sequences_label` - Gap a supplied set of sequences according to the `IMGT alignment <https://www.imgt.org/IMGTindex/AlleleAlignments.php>`_, warning of any missing conserved residues. It employs a gapped reference set, which must be provided. For each sequence to be gapped, the tool identifies the closest reference sequence, and uses it as a template.

:ref:`fix_macaque_gaps_label` - The IMGT alignment for macaque IG sequences has inserted codons relative to the alignment used for most other species. This can cause problems for downstream tools. This utility removes the additional codons in macaque IG, reverting to the standard alignment

Using custom databases with IgBlast
-----------------------------------

These tools are covered in detail in :ref:`the next section <custom_igblast_label>`.

:ref:`make_igblast_ndm_label` - Create a custom ndm file for IgBlast

:ref:`annotate_j_label` - Create a custom auxiliary data file for IgBlast

Convenience Tools
-----------------

:ref:`extract_refs_label` - Download reference sequences from IMGT for a specific species and locus

:ref:`identical_seqs_label` - Report cases where, in a FASTA file, the same sequence is listed more than once with different IDs

:ref:`rev_comp_label` - Reverse-complement a nucleotide sequence

:ref:`at_coords_label` - Find the sequence at specific co-ordinates within the sequence in a fasta file

:ref:`merge_fasta_label` - Merge sequences in two FASTA files, removing duplicates


simple_bio_seq API
------------------

This is a set of simple wrappers around commonly-used functions in Biopython which suit my use case and are used in the tools above. They support, for example,
one-line reading and writing of FASTA files and comma/tab-separated files, and manage sequences as strings in a dict, which I find simplifies code compared to the
Biopython functions, at the expense of flexibility which I rarely need. They are documented under :ref:`API <api_label>`.