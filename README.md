# receptor_utils

Some tools I find useful for working with Ig receptor sequences.

# Installation

```bash
pip install receptor-utils
```

The module requires [Biopython](https://biopython.org).

# Overview

Please refer to the files themselves for slightly more detailed documentation.

### simple_bio_seq 

Contains some convenience functions that are backed by BioPython but simplified for 
my use case. It uses the following approach to keep things simple
(at the expense of some flexibility/scalability):

- store sequences as strings, use dicts for collections
- convert sequences to upper case on input
- coerce iterators into lists for ease of debugging


```python
from receptor_utils import simple_bio_seq as simple
seqs = simple.read_fasta('seqfile.fasta')  # read sequences into a dict with names as keys
seq = simple.read_single_fasta('seqfile.fasta')  # reads the first or only sequence into a string
seq = simple.reverse_complement(seq)
```
See the file for other functions.

### novel_allele_name

Contains the function ```name_novel()```, which will generate a name for a 'previously undocumented'
allele, given its sequence. The name will consist of the name of the nearest allele in a 
reference set provided to the function, suffixed by the SNPs that differentiate it,
for example:

IGHV1-69*01_a29g_c113t

Numbering of V-sequences uses the IMGT alignment. The naming convention follows that used by
[Tigger](https://tigger.readthedocs.io/en/stable/) and 
[VDJbase](https://vdjbase.org).

### number_ighv

Contains various functions for working with V-sequences according to the IMGT numbering scheme.
The most useful is ```gap_sequence()``` which will gap the provided V-sequence by using the closest sequence in a reference
set as a template.

# Example scripts

These may be useful in their own right, but also show how to use some of the functions 
mentioned above. Once the package is installed, you should be able to run these at the command
line without the .py extension, for example type
```shell
$ extract_refs --help
```
for help

### rev_comp
Return the reverse-complement of the specified nucleotide sequence.

### name_allele
Return a 'tigger style' name for an allele sequence (reference allele name suffixed with SNPs)
given a reference set

### extract_refs
A script which uses ```simple_bio_seq``` to extract 
files for particular loci and species from an IMGT reference file. 

### identical_seqs
A script which uses ```simple_bio_seq``` to list identical sequences and sub-sequences in a fasta file.

### gap_inferred
A script which will gap a set of sequences listed in a FASTA file, using the closest sequences
discovered from a reference set. The script will do its best to warn of issues with the reference sequences and with 
the gapped sequences it provides: please use the warnings to check that things are ok.

Sequences to be gapped are assumed to be complete at the 5' end. If necessary they should be gapped with dots
at the 5' end, so that the first nucleotide is at the correct position in the full-length sequence (in exactly
the same way that IMGT puts dots at the start of a reference sequence that is incomplete at the 5' end)

### make_igblast_ndm
A script which uses a set of IMGT-gapped V-sequences to create the ndm file 
[required by IgBLAST for a custom organism](https://ncbi.github.io/igblast/cook/How-to-set-up.html)

### annotate_j
Given a set of J sequences, identify the correct frame and location of the CDR3 end, by searching for 
the GxG motif.











