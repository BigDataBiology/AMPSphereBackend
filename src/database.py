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

gmsc_metadata = pd.read_csv('data/tables/GMSCMetadata.tsv',
            sep='\t', index_col=0,
            # These do not get guessed correctly
            dtype={
                'geographic_location': str,
                'environment_material': str,
                },
        )
gmsc_metadata['geographic_location'].replace(np.nan, '', inplace=True)
gmsc_metadata['environment_material'].replace(np.nan, '', inplace=True)
gmsc_metadata['specI'].replace(np.nan, None, inplace=True)
gmsc_metadata['microbial_source_p'].replace(np.nan, None, inplace=True)
gmsc_metadata['microbial_source_c'].replace(np.nan, None, inplace=True)
gmsc_metadata['microbial_source_o'].replace(np.nan, None, inplace=True)
gmsc_metadata['microbial_source_f'].replace(np.nan, None, inplace=True)
gmsc_metadata['microbial_source_g'].replace(np.nan, None, inplace=True)
gmsc_metadata['microbial_source_s'].replace(np.nan, None, inplace=True)
number_genes_per_amp = gmsc_metadata.value_counts('AMP').to_dict()
amp2gmsc = gmsc_metadata[['AMP']].groupby('AMP').groups

amps = pd.read_csv(
        'data/tables/AMP.tsv',
        sep='\t',
        index_col=0,
        dtype={
            'Antifam': 'category',
            'RNAcode': 'category',
            'metaproteomes': 'category',
            'coordinates': 'category',
        },
        )
amps.rename(columns={'Coordinates': 'coordinates'}, inplace=True)
