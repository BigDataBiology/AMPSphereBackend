import polars as pl
import pandas as pd
import numpy as np
import sqlite3

gtdb_taxon_to_rank = \
        pd.read_table('data/tables/GTDBTaxonRank.tsv',
                      index_col=0, usecols=(1,2)
                  ).squeeze().to_dict()


selectedcoprediction = \
        pd.read_table('data/tables/AMP_coprediction_AMPSphere.tsv.xz',
            index_col=0
        ).to_dict('index')


db = sqlite3.connect('./data/ampsphere_main_db/AMPSphere_latest.sqlite', check_same_thread=False)

GMSC_METADATA_COLUMNS = [
        ('accession', str),
        ( 'gene_sequence', str),
        ( 'AMP', str),
        ( 'sample', str),
        ( 'specI', str),
        ( 'is_metagenomic', str),
        ( 'geographic_location', str),
        ( 'latitude', str),
        ( 'longitude', str),
        ( 'general_envo_name', str),
        ( 'environment_material', str),
        ( 'microbial_source_d', str),
        ( 'microbial_source_p', str),
        ( 'microbial_source_c', str),
        ( 'microbial_source_o', str),
        ( 'microbial_source_f', str),
        ( 'microbial_source_g', str),
        ( 'microbial_source_s', str),
        ]

def _float_or_None(x):
    return float(x) if x else None

def _make_gmsc_metadata_df(data):
    return pl.DataFrame(data, schema=GMSC_METADATA_COLUMNS).with_columns([
        (pl.col('is_metagenomic') == "True").alias('is_metagenomic'),
        (pl.col('latitude').map_elements(_float_or_None, return_dtype=pl.Float32)).alias('latitude'),
        (pl.col('longitude').map_elements(_float_or_None, return_dtype=pl.Float32)).alias('longitude'),
        ])

_number_genes_per_amp = {}
def number_genes_per_amp(accession):
    if accession not in _number_genes_per_amp:
        [(n,)] = db.execute('SELECT COUNT(*) FROM Metadata WHERE AMP = ?;', [accession]).fetchall()
        _number_genes_per_amp[accession] = n
    return _number_genes_per_amp[accession]


amps = pl.read_csv(
        'data/tables/AMP.tsv',
        separator='\t',
        dtypes={
            'Antifam': pl.Categorical,
            'RNAcode': pl.Categorical,
            'metaproteomes': pl.Categorical,
            'coordinates': pl.Categorical,
        },
        )
amps = amps.rename({'Coordinates': 'coordinates'})
amps = amps.with_columns([
    (pl.col('Antifam') == 'Passed').alias('Antifam'),
    (pl.col('metaproteomes') == 'Passed').alias('metaproteomes'),
    (pl.col('metatranscriptomes') == 'Passed').alias('metatranscriptomes'),
    (pl.col('coordinates') == 'Passed').alias('coordinates'),
    ])

