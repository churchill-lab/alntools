===============================
alntools
===============================


.. image:: https://img.shields.io/pypi/v/alntools.svg
    :target: https://pypi.python.org/pypi/alntools

.. image:: https://img.shields.io/travis/mattjvincent/alntools.svg
    :target: https://travis-ci.org/mattjvincent/alntools

.. image:: https://readthedocs.org/projects/alntools/badge/?version=latest
    :target: https://alntools.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://pyup.io/repos/github/mattjvincent/alntools/shield.svg
    :target: https://pyup.io/repos/github/mattjvincent/alntools/
    :alt: Updates


**alntools** processes next-generation sequencing read alignments into a sparse compressed incidence matrix (aka Equivalence Classes) and stores it in a pre-defined binary format for efficient downstream analyses and storage. It enables us to compare, contrast, or combine the results of different alignment strategies.

* Free software: Apache Software License 2.0
* Documentation: https://churchill-lab.github.io/alntools/


Features
--------

* ``split`` divides a large bam file into smaller ones
* ``bam2ec`` preprocesses a bam file into binary equivalence class (EC) format
* ``bam2emase`` preprocesses a bam file into EMASE format
* ``ec2emase`` converts binary EC file into EMASE format
* ``emase2ec`` converts EMASE format into binary EC format
* ``range`` finds effective lengths of target sequences from alignment data

We are planning to add more features for preprocessing different formats of NGS read alignments at the population level and summarizing useful information from them. Stay tuned!

