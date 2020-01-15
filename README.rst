.. image:: https://travis-ci.org/cgat-developers/cgat-apps.svg?branch=master
    :target: https://travis-ci.org/cgat-developers/cgat-apps

===========
Project 057
===========

This repository contains R scripts to accompany:
Ellender TJ, Avery SV, Mahfooz K, et al. Embryonic progenitor pools generate diversity in fine-scale excitatory cortical subnetworks. Nat Commun. 2019;10(1):5224. Published 2019 Nov 19. doi:10.1038/s41467-019-13206-1


Installation
============

1. Install https://github.com/cgat-developers/cgat-core and https://github.com/cgat-developers/cgat-apps
2. Download the project057 github repository
3. Run ```setup.py develop``` in the project057 directory

Usage
=====

The scripts require a counts table, generated e.g. by featurecounts as input.

Run the ``cgat --help`` command to see what scripts are available and how to use them.
For example, to run multiple normalisations on data obtained from featurecounts

   cgat-singlecell normalisation --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv

