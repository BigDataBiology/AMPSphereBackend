from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship
from src.database import Base


class AMP(Base):
    """

    """
    __tablename__ = "AMP"
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
    Antifam = Column(String)
    RNAcode = Column(String)
    metaproteomes = Column(String)
    metatranscriptomes = Column(String)
    coordinates = Column(String)


class GMSCMetadata(Base):
    """
    gmsc	amp	sample	source	specI	is_metagenomic	geographic_location	latitude	longitude	general envo name	environment_material
    """
    __tablename__ = "Metadata"
    GMSC_accession = Column(String, primary_key=True, index=True)
    gene_sequence = Column(String)
    AMP = Column(String, ForeignKey("AMP.accession"), index=True)
    sample = Column(String, index=True)
    # TODO can we unify all microbial sources in this table into a single rank? maybe species level?
    microbial_source = Column(String, index=True)
    specI = Column(String, index=True)
    is_metagenomic = Column(Boolean, index=True)
    geographic_location = Column(String, index=True)
    latitude = Column(Float)
    longitude = Column(Float)
    general_envo_name = Column(String, index=True)
    environment_material = Column(String)


class GTDBTaxonRank(Base):
    """

    """
    __tablename__ = "GTDB_taxon_rank"
    scientific_name = Column(String, primary_key=True, index=True)
    rank = Column(String, index=True)


class GTDBTaxonomy(Base):
    """
    TODO can we simplify this table? It's really a lot and maybe cannot finish a query soon.
    """
    __tablename__ = "GTDB_taxonomy_tree"
    Domain = Column(String, index=True)
    Kingdom = Column(String, index=True)
    Phylum = Column(String, index=True)
    Class = Column(String, index=True)
    Order = Column(String, index=True)
    Family = Column(String, index=True)
    Genus = Column(String, index=True)
    Species = Column(String, index=True)
    Subspecies = Column(String, index=True)


class Statistics(Base):
    """
    To be updated.
    """
    __tablename__ = "Statistics"
    gmsc = Column(Integer, primary_key=True)
    amp = Column(Integer)
    sample = Column(Integer)
    source = Column(Integer)
    specI = Column(Integer)
    is_metagenomic = Column(Integer)
    geographic_location = Column(Integer)
    latitude = Column(Integer)
    longitude = Column(Integer)
    general_envo_name = Column(Integer)
    environment_material = Column(Integer)
    family = Column(Integer)
