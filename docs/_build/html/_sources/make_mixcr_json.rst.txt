.. _make_mixcr_json:

make_mixcr_json
================

.. argparse::
   :filename: ../src/scripts/make_mixcr_json.py
   :func: get_parser
   :prog: make_mixcr_json

Notes
-----

This command will create a germline reference set in `Mixcr's JSON format <https://https://mixcr.com/mixcr/reference/ref-repseqio-json-format/?h=json>`_ 
using as input a reference set in `AIRR-C MIAIRR format <https://docs.airr-community.org/en/latest/datarep/miairr.html>`_.
If species name and/or taxon id are not provided, the utility will attempt to determine them from the input file.
