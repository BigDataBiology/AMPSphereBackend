#!/usr/bin/env bash


TMP_DIR=tmp
ZENODO_ID=4606582

mkdir $TMP_DIR

./scripts/download_files.sh $TMP_DIR $ZENODO_ID

./scripts/create_gene_source_hrenvo.py

./scripts/generate_features.py

./scripts/generate_main_db_tables.py

sqlite3 data/ampsphere_main_db/AMPSphere_latest.sqlite < ./scripts/import.sql

./script/generate_mmseqs_search_db.sh

./script/generate_hmm_profiles_db.sh

./script/generate_pre_computed_family_data.py

python -m pytest tests/testing.py