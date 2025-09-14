.. _create_alignment_label:

create_alignment
================

Create a formatted alignment display from gap-aligned IG/TR alleles. This tool takes FASTA files 
containing IMGT-gapped nucleotide sequences and produces a human-readable alignment showing nucleotide 
sequences with amino acid translations and silent/nonsilent mutations. The tool should be used
on a set of alleles from the same gene.

The tool supports V, D, and J sequence types with type-specific formatting:

- **V sequences**: Display both nucleotides and amino acids with CDR region delineation
- **D sequences**: Display only nucleotides (no amino acid translations)  
- **J sequences**: Display both nucleotides and amino acids without CDR regions (for translation, the alleles should be gap-aligned to a codon boundary at the 5' end)

For V sequences, CDR regions are marked based on configurable codon coordinates. 

Sequences are automatically sorted alphabetically by allele name for consistent output formatting.

.. argparse::
   :filename: ../src/scripts/create_alignment.py
   :func: get_parser
   :prog: create_alignment

CDR delineation
---------------

For V sequences, the CDR regions are marked using codon coordinates (because the input sequences are IMGT-aligned, the CDR positions occupy specific coordinates). 
By default, the tool uses the following IMGT codon positions, which are correct for IMGT-aligned human allele sequences:
- CDR1 start: 27
- CDR1 end: 38
- CDR2 start: 56
- CDR2 end: 65
- CDR3 start: 105

Codon boundaries for other species may differ, and the user can specify alternative coordinates using the ``--v_coords`` option. As an example, the
default coordinates would be specified as follows::

    --v_coords "27,38,56,65,105"

Examples
--------

Create an alignment of V sequences with default CDR coordinates::

    create_alignment sequences.fasta V output.txt

Create an alignment of D sequences (nucleotides only)::

    create_alignment d_sequences.fasta D d_alignment.txt

Create an alignment with custom CDR coordinates and shorter line wrapping::

    create_alignment sequences.fasta V output.txt --v_coords "25,30,54,60,100" --codon_wrap 15

Filter sequences containing a specific pattern::

    create_alignment sequences.fasta V output.txt --filter "IGHV1-18"