from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship
from src.database import Base


class AMP(Base):
    """
accession       sequence        family  length  molecular_weight        isoelectric_point       charge  aromaticity     instability_index       gravy   Antifam RNAcode metaproteomes   metatranscriptomes
      coordinates
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
    def __str__(self):
        return f'AMP(accession={self.accession}, family={self.family})'
    __repr__ = __str__


class GMSCMetadata(Base):
    """
    accession       gene_sequence   AMP     sample  source  specI   is_metagenomic  geographic_location
    latitude        longitude       general envo name       environment_material    microbial_source_d
    microbial_source_p      microbial_source_c      microbial_source_o      microbial_source_f      microbial_source_g      microbial_source_s
    """
    __tablename__ = "Metadata"
    GMSC_accession = Column(String, primary_key=True, index=True)
    gene_sequence = Column(String)
    AMP = Column(String, ForeignKey("AMP.accession"), index=True)
    sample = Column(String, index=True)
    specI = Column(String, index=True)
    is_metagenomic = Column(Boolean, index=True)
    geographic_location = Column(String, index=True)
    latitude = Column(Float)
    longitude = Column(Float)
    general_envo_name = Column(String, index=True)
    environment_material = Column(String)
    microbial_source_d = Column(String, index=True)  # kingdom
    microbial_source_p = Column(String, index=True)  # phylum
    microbial_source_c = Column(String, index=True)  # class
    microbial_source_o = Column(String, index=True)  # order
    microbial_source_f = Column(String, index=True)  # family
    microbial_source_g = Column(String, index=True)  # genus
    microbial_source_s = Column(String, index=True)  # species


class GTDBTaxonRank(Base):
    """

    """
    __tablename__ = "GTDBTaxonRank"
    record_idx = Column(Integer, primary_key=True)
    gtdb_taxon = Column(String, index=True)
    microbial_source_rank = Column(String, index=True)


class Statistics(Base):
    """
    To be updated.
    """
    __tablename__ = "Statistics"
    gmsc = Column(Integer, primary_key=True)
    amp = Column(Integer)
    sample = Column(Integer)
    source_d = Column(Integer)  # kingdom
    source_p = Column(Integer)  # phylum
    source_c = Column(Integer)  # class
    source_o = Column(Integer)  # order
    source_f = Column(Integer)  # family
    source_g = Column(Integer)  # genus
    source_s = Column(Integer)  # species
    geographic_location = Column(Integer)
    general_envo_name = Column(Integer)
    environment_material = Column(Integer)
    sphere = Column(Integer)
    family = Column(Integer)
