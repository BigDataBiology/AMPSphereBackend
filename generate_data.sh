#!/usr/bin/env bash


TMP_DIR=tmp
VERSION_CODE='2022-03'
mkdir -p $TMP_DIR

set -ev

# Download necessary files
# Download GTDB files
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_taxonomy_r95.tsv.gz -P $TMP_DIR --no-check-certificate
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_metadata_r95.tar.gz -P $TMP_DIR --no-check-certificate
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_metadata_r95.tar.gz -P $TMP_DIR --no-check-certificate
wget -c https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_taxonomy_r95.tsv.gz -P $TMP_DIR --no-check-certificate

# Generate necessary tables in order to create a sqlite database.
python ./scripts/generate_main_db_tables.py ${TMP_DIR} \
    --metadata data_folder/gmsc_amp_genes_envohr_source.tsv.gz \
    --faa data_folder/AMPSphere_v.${VERSION_CODE}.faa.gz \
    --fna data_folder/AMPSphere_v.${VERSION_CODE}.fna.xz \
    --features data_folder/AMPSphere_v.${VERSION_CODE}.features_for_web.tsv.gz \
    --quality data_folder/quality_assessment.tsv.gz \
    --gtdb-files ${TMP_DIR}/bac120_metadata_r95.tar.gz/bac120_metadata_r95.tsv ${TMP_DIR}/ar122_metadata_r95.tar.gz/ar122_metadata_r95.tsv

# Move generated tables
mkdir -p data/tables
mv ${TMP_DIR}/AMP.tsv data/tables/AMP.tsv
mv ${TMP_DIR}/GMSCMetadata.tsv data/tables/GMSCMetadata.tsv
mv ${TMP_DIR}/GTDBTaxonRank.tsv data/tables/GTDBTaxonRank.tsv
mv ${TMP_DIR}/Statistics.tsv data/tables/Statistics.tsv

# Generate sqlite3 database
DB_DIR='data/ampsphere_main_db'
mv ${DB_DIR}/AMPSphere_latest.sqlite ${DB_DIR}/AMPSphere_latest.sqlite.backup
python -c "from src.database import engine; from src import models; models.Base.metadata.create_all(bind=engine)"
sqlite3 ${DB_DIR}/AMPSphere_latest.sqlite < ./scripts/import.sql

# Copy precomputed family data
mkdir -p data/pre_computed/families
cp -r data/original_data/AMPSphere_generation_v.${VERSION_CODE}/analysis/families/* data/pre_computed/families

# Generate mmseqs database and index
# TODO Tested
DB_DIR='data/mmseqs_db'
mkdir -p $DB_DIR
mmseqs createdb data_folder/AMPSphere_v.${VERSION_CODE}.faa.gz \
    ${DB_DIR}/AMPSphere_latest.mmseqsdb
mmseqs createindex ${DB_DIR}/AMPSphere_latest.mmseqsdb ${DB_DIR}/tmp

# Generate hmmprofile database
DB_DIR='data/hmmprofile_db'
rm -rf $DB_DIR
mkdir -p $DB_DIR
cat data/pre_computed/families/hmm/*.hmm > data/hmmprofile_db/AMPSphere_latest.hmm
hmmpress data/hmmprofile_db/AMPSphere_latest.hmm
