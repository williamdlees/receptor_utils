.. _airrc_sets_with_MiXCR:

Using AIRR Community Reference Sets with MiXCR
==============================================

The utility :ref:`make_mixcr_json` will convert AIRR-C germline reference sets in MIAIRR format to MiXCR's JSON format. In addition,
the :ref:`download_germline_set` utility can be used to download germline reference sets from the AIRR-C OGRDB database in MiXCR format. 
The germline set in MiXCR format can then be used to annotate sequences with MiXCR, using `mixcr analyze` as outlined in the MiXCR
documentation on the `Custom Reference Library page <https://mixcr.com/mixcr/guides/create-custom-library/?h=json#de-novo-libraries>`_.
There is no need to use `mixcr buildlibrary` as the file is already in the correct format for use with MiXCR.

Example Usage
=============

For simplicity in this example, we will use the MiXCR Docker/Singularity container.  
If you are not familiar with the container, you can read about it in the `MiXCR documentation <https://mixcr.com/mixcr/getting-started/docker/?h=container>`_. Alternatively, you can install MiXCR on your local machine,
using the instructions provided on the `MiXCR website <https://mixcr.com>`_. In this case just skip the 'log in to Docker' step in the walkthrough.

We will annotate a sample set of sequences provided on the `Immcantation website <http://clip.med.yale.edu/immcantation/examples/AIRR_Example.tar.gz>`_. 
The sequences have been preprocessed and quality-filtered from Illumina paired-end reads and are present in a FASTA file called ``HD13M.fasta``. The process required to create such a file from sequencing reads, will depend on the sequencing protocol. 
You can consult the Mixcr documentation, or other sources, to determine a suitable approach for you sequencing data.

Prerequisites
-------------
Before you start, you will need to have the following installed on your machine:

    * Docker
    * Python 3.9 or above
    * The receptor-utils package (see :ref:`introduction_label` for installation instructions)
    * The file ``HD13M.fasta``, extracted from the tarball which you can download `here <http://clip.med.yale.edu/immcantation/examples/AIRR_Example.tar.gz>`_.
  

The steps we shall follow are as follows:

#. Download the AIRR-C germline set for human IGHV genes in Mixcr format.
#. Log in to the MiXCR container so that we can use its installed tools.
#. Annotate the sample sequences using MiXCR.

1. Download the germline set for human IGHV genes
-------------------------------------------------

In a suitable directory to use for this test, use the :ref:`download_germline_set` utility to download the germline set for human IGHV genes:

..  code-block:: none

    $ download_germline_set "Homo sapiens" IGH -n IGH_VDJ -f MIXCR
    https://ogrdb.airr-community.org/api_v2/germline/species
    Homo sapiens: 9606
    9606.IGH_VDJ
    MIXCR JSON file saved to Homo_sapiens_IGH_VDJ_IGH_mixcr.json
    
Finally, extract the file ``HD13M.fasta`` from the downloaded tarball and copy it to the directory.

1. Log in to the MiXCR container so that we can use its installed tools
-----------------------------------------------------------------------

Log in to the container, mounting the current local directory as /work.

From Linux using Docker:

.. code-block:: none

    docker run -it -v $(pwd):/work ghcr.io/milaboratory/mixcr/mixcr:latest bash


From Windows using Docker:

.. code-block:: none

    docker run -it -v %cd%:/work ghcr.io/milaboratory/mixcr/mixcr:latest bash


Once in the container, cd to /work and check that the reference set files are present:

.. code-block:: none

    bash-4.2# cd /work
    bash-4.2# ls
    HD13M.fasta Homo_sapiens_IGH_VDJ_IGH_mixcr.json
    bash-4.2#

3. Set the MiXCR license to match your key:
-------------------------------------------

Set the MiXCR license to match your key:

.. code-block:: none

    bash-4.2# MI_LICENSE="...your licence key here..."
    bash-4.2# export MI_LICENSE


4. Annotate the sample sequences using MiXCR
--------------------------------------------

Annotate the sequences in ``HD13M.fasta`` with MiXCR. As the sequences are in plain FASTA format, we will use the 'generic pacbio' template.

.. code-block:: none

	bash-4.2# mixcr analyze generic-pacbio -s "Homo sapiens" \
        --library Homo_sapiens_IGH_VDJ_IGH_mixcr \
        --assemble-clonotypes-by FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4 \
        HD13M.fasta \
        HD13M



Once MiXCR has run, you can view ``HD13M.clones_IGH.tsv`` to see the resulting annotations.

