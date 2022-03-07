.mode csv
.separator \t
.import --skip 1 tmp/AMP.tsv AMP
.import --skip 1 tmp/GMSCMetadata.tsv Metadata
.import --skip 1 tmp/GTDBTaxonRank.tsv GTDBTaxonRank
.import --skip 1 tmp/Statistics.tsv Statistics
