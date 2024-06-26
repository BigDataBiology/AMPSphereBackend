import polars as pl
import pandas as pd
import numpy as np

gtdb_taxon_to_rank = \
        pd.read_table('data/tables/GTDBTaxonRank.tsv',
                      index_col=0, usecols=(1,2)
                  ).squeeze().to_dict()

coprediction = \
        pd.read_table('data/tables/AMP_coprediction_AMPSphere.tsv.xz',
            index_col=0
        ).to_dict('index')

gmsc_metadata = pl.read_csv(
        'data/tables/GMSCMetadata.tsv',
        separator='\t',
        dtypes={
            'latitude': pl.Float32,
            'longitude': pl.Float32,
            },
        )

counts = gmsc_metadata.group_by('AMP').len()
number_genes_per_amp = dict(zip(counts['AMP'], counts['len']))
del counts

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

