##receptor_utils

Some tools I find useful for working with Ig receptor sequences.

## Installation

git clone https://github.com/williamdlees/receptor_utils
pip install receptor_utils

The module requires Biopython to be installed.

(will be on PyPi soon)

## Overview

# simple_bio_seq 

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

# novel_allele_name

Contains the function name_novel(), which will generate a name for a 'previously undocumented'
allele, given its sequence. The name will consist of the name of the nearest allele in a 
reference set provided to the function, suffixed by the SNPs that differentiate it,
for example:

IGHV1-69*01_a29g_c113t

Numbering of V-sequences uses the IMGT alignment. The naming convention follows that used by
[Tigger](https://tigger.readthedocs.io/en/stable/) and 
[VDJbase](https://vdjbase.org).

# number_ighv

Contains various functions for working with  V-sequences according to the IMGT numbering scheme.
The most useful is ```gap_sequence()``` which will gap the provided V-sequence by using the closest sequence in a reference
set as a template.

# extract_refs.py
A script (that  can be run from the cloned repo) which uses simple_bi_seq to extract 
files for particular loci and species from an IMGT reference file. 

# extract_refs.py
A script (that can be run from the cloned repo) which uses simple_bi_seq to extract 
list identical sequences and sub-sequences in a fasta file.














