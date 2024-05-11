from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
from src.utils import cfg
import pandas as pd

SQLALCHEMY_DATABASE_URL = "sqlite:///" + cfg['ampsphere_main_db']

engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()

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
number_genes_per_amp = gmsc_metadata.value_counts('AMP').to_dict()

amps = pd.read_csv('data/tables/AMP.tsv', sep='\t', index_col=0)
amps.rename(columns={'Coordinates': 'coordinates'}, inplace=True)
