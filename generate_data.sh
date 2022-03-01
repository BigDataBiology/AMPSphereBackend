#!/usr/bin/env bash


TMP_DIR=tmp
ZENODO_ID=4606582
VERSION_CODE='2021-03'
mkdir $TMP_DIR


# Download necessary files
# Quality table should be directly included in the extra input files as later we'll need to include the file in Zenodo.
# TODO Tested
./scripts/download_files.sh ${TMP_DIR} ${ZENODO_ID}

# Create gmsc_amp_genes_envohr_source.tsv
# Help needed: Can you prepare a script to generate the file? It uses a lot of scripts and multiple data files, it's hard to understand these scripts again.
./scripts/generate_gmsc_amp_genes_envohr_source.sh

# Create feature tables for all the AMPs.
# TODO Tested
./scripts/generate_features.py

# Generate necessary tables in order to create a sqlite database.
DB_DIR='data/ampsphere_main_db'
./scripts/generate_main_db_tables.py

# Generate sqlite3 database
sqlite3 ${DB_DIR}/AMPSphere_latest.sqlite < ./scripts/import.sql

# Precompute family data
# Help needed: Please give me the script.
./script/pre_compute_family_data.py

# Generate mmseqs database and index
# TODO Tested
DB_DIR='data/mmseqs_db'
./script/generate_mmseqs_search_db.sh ${TMP_DIR}/AMPSphere_v.${VERSION_CODE}.faa.gz ${DB_DIR}/AMPSphere_latest.sqlite ${DB_DIR}/tmp

# Generate hmmprofile database
# Help needed: Please give me the script, but maybe the file were already generated from `pre_compute_family_data.py`?
DB_DIR='data/hmmprofile_db'
./script/generate_hmmprofiles_db.sh

# Do unit test
python -m pytest tests/testing.py

