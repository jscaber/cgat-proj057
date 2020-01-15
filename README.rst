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
2. Install https://github.com/AllenInstitute/scrattch.hicat
3. Download the project057 github repository
4. Run ```setup.py develop``` in the project057 directory

Usage
=====

The scripts require a ```.tsv``` counts table, generated e.g. by featurecounts as input.
They also require an annotation table for the columns.

Run the ``cgat-singlecell --help`` command to see what scripts are available and how to use them.

For example, to run filter data obtained from featurecounts
```cgat-singlecell filter --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv```

To run the normalisation script and map the data onto a reference dataset (e.g. Allen) use:
```cgat-singlecell normalisation --rds-filename sce_filtered_hicat.rds --ERCC ERCC.tsv --allen-design mouse_VISp_2018-06-14_samples-columns.csv --allen-datamatrix mouse_VISp_2018-06-14_exon-matrix.csv --allen-rowdata mouse_VISp_2018-06-14_genes-rows.csv --allen-filter "L4 IT" --norm scran --colours red,green --allen-colours black --perplexity 5 > pipeline.log```

Analysis Sequence
=================

For the publication the following sequence was used on each experimental plate separately.

1. Filtering of datase using ```cgat-singlecell filtering```
2. Run scrattch.hicat using wrapper script ```cgat-singlecell hicat```
3. Normalisation and mapping to reference dataset using ```cgat-singlecell normalisation```
4. Differential expression analysis using ```cgat-singlecell sc-diffexpression```


