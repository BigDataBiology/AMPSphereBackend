.mode csv
.separator \t
.import --skip 1 data/tables/AMP.tsv AMP
.import --skip 1 data/tables/GMSCMetadata.tsv Metadata
.import --skip 1 data/tables/GTDBTaxonRank.tsv GTDBTaxonRank
.import --skip 1 data/tables/Statistics.tsv Statistics
