import shutil
import gzip
import pandas as pd
import argparse
from Bio import SeqIO
import datatable as dt
import pathlib


parser = argparse.ArgumentParser()
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
metadata = dt.fread(args.metadata, sep='\t').to_pandas()
features = dt.fread(args.features, sep='\t').to_pandas()
with gzip.open(args.faa, 'rt') as f:
    faa = SeqIO.parse(f, 'fasta')
    AMP_table = pd.DataFrame([[r.id, str(r.seq), r.description.split(' | ')[1]] for r in faa.records],
                             columns=['accession', 'sequence', 'family'])
with gzip.open(args.fna, 'rt') as f:
    fna = SeqIO.parse(f, 'fasta')
    gmsc_table = pd.DataFrame([[r.id, str(r.seq), r.description.split(' ')[1]] for r in fna.records],
                              columns=['accession', 'gene_sequence', 'AMP'])
quality = pd.read_csv('tmp/quality_assessment.tsv.xz', sep='\t')
taxonomy = pd.concat((pd.read_csv(file, sep='\t') for file in args.gtdb_files), axis=0)
# TODO generate hierarchical taxonomy tables based on selected field, but which field?

print('done')


output_dir = pathlib.Path(args.outdir)
if not output_dir.is_dir():
    output_dir.mkdir()


print('Generating tables...', end=' ')
"""
    accession = Column(String, primary_key=True, index=True)
    sequence = Column(String)
    family = Column(String, index=True)
    length = Column(Integer, index=True)
    molecular_weight = Column(Float, index=True)
    isoelectric_point = Column(Float, index=True)
    charge = Column(Float, index=True)
    aromaticity = Column(Float)
    instability_index = Column(Float)
    gravy = Column(Float)
"""
AMP_cols = ['accession', 'sequence', 'family', 'length', 'molecular_weight', 'isoelectric_point', 'charge', 'aromaticity', 'instability_index', 'gravy']

tables = dict(
    AMP=pd.merge(AMP_table.drop(columns='family'), features, left_on=['accession'], right_on=['id'])[AMP_cols].
                     merge(quality, left_on='accession', right_on='AMP'),
    GMSC=gmsc_table,
    Statistics=pd.DataFrame({**metadata.nunique(dropna=True).to_frame().T.to_dict(),
                                **AMP_table.drop(columns=['sequence', 'accession']).nunique(dropna=True).to_frame().T.to_dict()}))

# SPHEREs with num_amps < 8 should not be treated as families.
tables['Statistics']['family'] = sum(tables['AMP'].family.value_counts() >= 8)
print('done')

shutil.copy(metadata_file, pathlib.Path(output_dir.joinpath('Metadata.tsv')))
for name, table in tables.items():
    outfile = output_dir.joinpath(name + '.tsv')
    print('Writing to {}...'.format(outfile), end=' ')
    table.to_csv(outfile, sep='\t', index=False, header=False)
    print('done')
