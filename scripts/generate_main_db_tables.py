import shutil
import gzip
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import pathlib
import xz
from collections import defaultdict


parser = argparse.ArgumentParser(description='''
Under project root dir, run
python scripts/generate_main_db_tables.py tmp \
--metadata tmp/gmsc_amp_genes_envohr_source.tsv \
--faa tmp/AMPSphere_v.2021-03.faa.gz \
--fna tmp/AMPSphere_v.2021-03.fna.xz \
--features tmp/features_plot_ampsphere.tsv \
--quality tmp/quality_assessment.tsv.xz \
--gtdb-files tmp/bac120_metadata_r95.tar.gz/bac120_metadata_r95.tsv tmp/ar122_metadata_r95.tar.gz/ar122_metadata_r95.tsv
''')
parser.add_argument('--metadata', type=str)
parser.add_argument('--faa', type=str)
parser.add_argument('--fna', type=str)
parser.add_argument('--features', type=str)
parser.add_argument('--quality', type=str)
parser.add_argument('--gtdb-files', type=str, nargs='+')
parser.add_argument('outdir', type=str, default='tables')
args = parser.parse_args()


print('Loading input data...', end=' ')
metadata_file = pathlib.Path(args.metadata)
metadata = pd.read_table(args.metadata, sep='\t', na_values=['NA'])
features = pd.read_table(args.features, sep='\t')
with gzip.open(args.faa, 'rt') as f:
    faa = SeqIO.parse(f, 'fasta')
    AMP_cols = ['accession', 'sequence', 'family']
    AMP_table = pd.DataFrame([[r.id, str(r.seq), r.description.split(' | ')[1]] for r in faa.records], columns=AMP_cols)
with xz.open(args.fna, 'rt') as f:
    fna = SeqIO.parse(f, 'fasta')
    GMSC_cols = ['accession', 'gene_sequence', 'AMP']
    gmsc_table = pd.DataFrame([[r.id, str(r.seq), r.description.split(' | ')[1]] for r in fna.records], columns=GMSC_cols)
quality = pd.read_table(args.quality, sep='\t')
taxa_files = [pd.read_table(file, sep='\t') for file in args.gtdb_files]
taxonomy = pd.concat(taxa_files).reset_index()
# Generate progenomes to gtdb conversion dict.
ncbi_id_to_gtdb_tax = defaultdict(str)
ncbi_id_to_gtdb_tax.update(dict(zip(taxonomy.ncbi_taxid.tolist(), taxonomy.gtdb_taxonomy.apply(lambda x: x.split(';')[-1].split('__')[-1]))))
taxonomy = taxonomy['gtdb_taxonomy']
# Generate hierarchical taxonomy tables based on selected field, but which field?
taxonomy = taxonomy.str.split(';', expand=True).applymap(lambda x: x.split('__')[1])
taxonomy_cols = ('microbial_source_' + pd.Series(list('dpcofgs'))).tolist()
column_rename = dict(zip(range(7), taxonomy_cols))
taxonomy = taxonomy.rename(columns=column_rename)
taxonomy = pd.DataFrame({tax: row[:col_idx + 1] for row_idx, row in taxonomy.iterrows() for col_idx, tax in enumerate(row.tolist())})
taxonomy = taxonomy.T[taxonomy_cols].reset_index()
print(taxonomy)


output_dir = pathlib.Path(args.outdir)
if not output_dir.is_dir():
    output_dir.mkdir()


print('Generating tables...', end=' ', flush=True)
print('\t AMP table...', end=' ', flush=True)

AMP_cols += ['length', 'molecular_weight', 'isoelectric_point', 'charge', 'aromaticity', 'instability_index', 'gravy']
AMP_table = AMP_table.drop(columns='family').merge(features, left_on=['accession'], right_on=['id'])[AMP_cols]
AMP_cols += ['Antifam', 'RNAcode', 'metaproteomes', 'metatranscriptomes', 'Coordinates']
AMP_table = AMP_table.merge(quality, left_on='accession', right_on='AMP')[AMP_cols]
AMP_table.to_csv(output_dir.joinpath('AMP.tsv'), sep='\t', index=False)
del AMP_table, features, quality

print('\t GMSC table...', end=' ', flush=True)
metadata['source'] = metadata.apply(lambda row: row['source'] if row['is_metagenomic']
                else ncbi_id_to_gtdb_tax[int(row['sample'].split('.')[0])], axis=1)
GMSC_cols += metadata.columns.tolist()
gmsc_table = gmsc_table.merge(metadata, left_on='accession', right_on='gmsc', how='outer')[GMSC_cols]
del metadata
GMSC_cols += taxonomy_cols
gmsc_table = gmsc_table.merge(taxonomy, left_on='source', right_on='index', how='left', copy=False)[GMSC_cols]
gmsc_table = gmsc_table.drop(columns=['gmsc', 'amp', 'source'])
gmsc_table.loc[gmsc_table.is_metagenomic == False, 'general_envo_name'] = 'Progenomes'
gmsc_table.to_csv(output_dir.joinpath('GMSCMetadata.tsv'), sep='\t', index=False)
del gmsc_table

print('\t GTDBTaxonRank table...', end=' ', flush=True)
taxonomy = taxonomy.set_index('index')
rank = pd.DataFrame([{'gtdb_taxon': taxon, 'microbial_source_rank': col} for col in taxonomy.columns for taxon in taxonomy[col] if taxon])
rank.drop_duplicates().dropna().reset_index().to_csv(output_dir.joinpath('GTDBTaxonRank.tsv'), sep='\t', index=False)

