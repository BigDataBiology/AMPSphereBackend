#!/usr/bin/env bash

# $1 ZENODO_ID
# $2 TMP_DIR

# Download the Zenodo repository
zenodo_get -r $1 -o $2

# Download the extra input files provided by Celio.
scp -r ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/resource_table/* $2
wget -c http://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab -O ${2}/pgenomes_samples.tsv

# Download GTDB files.
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_taxonomy_r95.tsv.gz -P $2 --no-check-certificate
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_metadata_r95.tar.gz -P $2 --no-check-certificate
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_metadata_r95.tar.gz -P $2 --no-check-certificate
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_taxonomy_r95.tsv.gz -P $2 --no-check-certificate
