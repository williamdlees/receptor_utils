# https://ao.ms/how-to-make-a-python-script-pip-installable/

from setuptools import setup

setup(
    name="receptor_utils",
    version="1.0.0",
    scripts=["extract_refs.py",
             "identical_seqs.py",
             "novel_allele_name.py",
             "number_ighv.py",
             "simple_bio_seq.py"]
)
