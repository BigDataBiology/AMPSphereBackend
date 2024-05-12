from typing import List, Optional, Dict
from pydantic import field_validator, ConfigDict, BaseModel


# Define JSON objects to be returned to frontend here.

# Plot objects ------------------------------------------------------
class LinePlotData(BaseModel):
    type: str
    x: List[str]
    y: List[float]
    c: List[str]
    model_config = ConfigDict(from_attributes=True)


class FeatureGraphPoints(BaseModel):
    transfer_energy: LinePlotData
    hydrophobicity_parker: LinePlotData
    surface_accessibility: LinePlotData
    flexibility: LinePlotData


class BarPlotData(BaseModel):
    type: str = 'bar plot'
    labels: List[str]
    values: List[float]
    model_config = ConfigDict(from_attributes=True)


class BubbleMapData(BaseModel):
    type: str = 'bubble map'
    lat: List[float]
    lon: List[float]
    size: List[float]
    colors: List[str] = ['']
    model_config = ConfigDict(from_attributes=True)


class Distributions(BaseModel):
    geo: BubbleMapData
    habitat: BarPlotData
    microbial_source: Optional[BarPlotData] = None
    model_config = ConfigDict(from_attributes=True)


class AMPFeatures(BaseModel):
    MW: float
    Length: float
    Molar_extinction: Dict[str, float]
    Aromaticity: float
    GRAVY: float
    Instability_index: float
    Isoelectric_point: float
    Charge_at_pH_7: float
    Secondary_structure: Dict[str, float]
    graph_points: Optional[FeatureGraphPoints] = None
    model_config = ConfigDict(from_attributes=True)


# Object for AMP_card page ------------------------------------------------
class Metadata(BaseModel):
    AMP: str
    GMSC_accession: str
    gene_sequence: str
    sample: str
    specI: Optional[str] = None
    is_metagenomic: bool
    geographic_location: str
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    general_envo_name: str
    environment_material: str
    microbial_source_d: Optional[str] = None
    microbial_source_p: Optional[str] = None
    microbial_source_c: Optional[str] = None
    microbial_source_o: Optional[str] = None
    microbial_source_f: Optional[str] = None
    microbial_source_g: Optional[str] = None
    microbial_source_s: Optional[str] = None
    model_config = ConfigDict(from_attributes=True)

    @field_validator(
            'specI',
            'latitude',
            'longitude',
            'microbial_source_d',
            'microbial_source_p',
            'microbial_source_c',
            'microbial_source_o',
            'microbial_source_f',
            'microbial_source_g',
            'microbial_source_s',
            mode="before")
    @classmethod
    def empty_string_to_null(cls, val):
        return None if val == "" else val


class PageInfo(BaseModel):
    currentPage: int
    pageSize: int
    totalPage: int
    totalItem: int
    model_config = ConfigDict(from_attributes=True)


class PagedMetadata(BaseModel):
    info: PageInfo
    data: List[Metadata]
    model_config = ConfigDict(from_attributes=True)


class AMPQuality(BaseModel):
    Antifam: str
    RNAcode: str
    metaproteomes: str
    coordinates: str
    score: float
    badge: str
    model_config = ConfigDict(from_attributes=True)


class AMP(BaseModel):
    accession: str
    sequence: str
    family: str
    length: int
    molecular_weight: float
    isoelectric_point: float
    charge: float
    aromaticity: float
    instability_index: float
    gravy: float
    Antifam: str
    RNAcode: str
    metaproteomes: str
    metatranscriptomes: str
    coordinates: str
    num_genes : Optional[int] = None
    secondary_structure: Optional[Dict[str, float]] = None
    metadata: Optional[PagedMetadata] = None
    model_config = ConfigDict(from_attributes=True)


class PagedAMPs(BaseModel):
    info: PageInfo
    data: List[AMP]
    model_config = ConfigDict(from_attributes=True)


# Object for Family page ------------------------------------------------
class FamilyDownloads(BaseModel):
    alignment: str
    sequences: str
    hmm_profile: str
    tree_figure: str
    tree_nwk: str
    model_config = ConfigDict(from_attributes=True)




class Family(BaseModel):
    accession: str
    consensus_sequence: str
    num_amps: int
    downloads: FamilyDownloads
    associated_amps: List[str]
    feature_statistics: Dict[str, AMPFeatures]
    distributions: Distributions
    model_config = ConfigDict(from_attributes=True)


class PagedFamilies(BaseModel):
    info: PageInfo
    data: List[Family]
    model_config = ConfigDict(from_attributes=True)


# Object for Download page ------------------------------------------------

class mmSeqsSearchResult(BaseModel):
    query_identifier: str
    target_identifier: str
    sequence_identity: float
    alignment_length: int
    number_mismatches: int
    number_gap_openings: int
    domain_start_position_query: int
    domain_end_position_query: int
    domain_start_position_target: int
    domain_end_position_target: int
    E_value: float
    bit_score: int
    seq_query: str
    seq_target: str
    alignment_strings: Optional[List[str]] = None
    family: Optional[str] = None
    model_config = ConfigDict(from_attributes=True)


class HMMERSearchResult(BaseModel):
    query_accession: str
    query_length: int
    query_name: str
    target_accession: str
    target_length: int
    target_name: str
    E_value: float
    acc: float
    bias: float
    c_Evalue: float
    i_Evalue: float
    num_domain: int
    domain_index: int
    score: float
    from_ali: int
    from_env: int
    from_hmm: int
    to_ali: int
    to_env: int
    to_hmm: int
    description_of_target: str
    model_config = ConfigDict(from_attributes=True)


class Statistics(BaseModel):
    num_genes: int
    num_amps: int
    num_families: int
    num_habitats: int
    num_genomes: int
    num_metagenomes: int
    model_config = ConfigDict(from_attributes=True)


class Filters(BaseModel):
    family: List[str]
    habitat: List[str]
    sample_genome: List[str]
    microbial_source: List[str]
    model_config = ConfigDict(from_attributes=True)
