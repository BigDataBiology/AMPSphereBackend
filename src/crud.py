import math
import pathlib
import types

import pandas as pd
import sqlalchemy.exc
from sqlalchemy.orm import Session, Query
from sqlalchemy import or_, and_, not_
from sqlalchemy import distinct, func, select
from sqlalchemy.sql.expression import func as sql_func
from fastapi import HTTPException
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
from Bio.pairwise2 import format_alignment


from src import models, utils, database

def filter_by_criteria(query, db, **criteria):
    """
    Possible criteria:
        exp_evidence: str = None,
        antifam: str = None,
        RNAcode: str = None,
        coordinates: str = None,
        family: str = None,
        habitat: str = None,
        sample: str = None,
        microbial_source: str = None,
        pep_length_interval: str = None,
        mw_interval: str = None,
        pI_interval: str = None,
        charge_interval: str = None,
    """
    cols_mapper = dict(   # Mapping from query filter names to table column names
        GMSCMetadata=dict(
            habitat='general_envo_name',
            microbial_source='microbial_source',
            sample_genome='sample'),
        AMP=dict(
            antifam='Antifam',
            RNAcode='RNAcode',
            coordinates='coordinates',
            pep_length_interval='length',
            mw_interval='molecular_weight',
            pI_interval='isoelectric_point',
            charge_interval='charge',
            family='family'))

    criteria = {key: value for key, value in criteria.items() if value}

    subquery = None
    for filter, value in criteria.items():
        if filter in {'habitat', 'sample_genome'}:
            if subquery is None:
                subquery = db.query(models.GMSCMetadata.AMP)
            subquery = subquery.where(getattr(models.GMSCMetadata, cols_mapper['GMSCMetadata'][filter]) == value)
        elif filter == 'microbial_source':
            if value not in database.gtdb_taxon_to_rank:
                raise HTTPException(status_code=400, detail='wrong taxonomy name provided.')
            rank = database.gtdb_taxon_to_rank[value]
            if subquery is None:
                subquery = db.query(models.GMSCMetadata.AMP)
            subquery = subquery.where(getattr(models.GMSCMetadata, rank) == value)
        elif filter == 'exp_evidence' and value == 'Passed':
            query = query.filter(or_(models.AMP.metaproteomes == 'Passed', models.AMP.metatranscriptomes == 'Passed'))
        elif filter == 'exp_evidence' and value == 'Failed':
            query = query.filter(and_(models.AMP.metaproteomes == 'Failed', models.AMP.metatranscriptomes == 'Failed'))
        elif filter in {'family', 'antifam', 'RNAcode', 'coordinates'}:
            query = query.filter(getattr(models.AMP, cols_mapper['AMP'][filter]) == value)
        elif filter in {'pep_length_interval', 'mw_interval', 'pI_interval', 'charge_interval'}:
            min_max: [str] = value.split(',')
            col_values = getattr(models.AMP, cols_mapper['AMP'][filter])
            query = query.filter(and_(float(min_max[0]) <= col_values, col_values <= float(min_max[1])))
    if subquery:
        query = query.where(models.AMP.accession.in_(subquery))
    return query



def get_amps(db: Session, page: int = 0, page_size: int = 20, **kwargs):
    query = db.query(models.AMP)
    query = filter_by_criteria(query=query, db=db, **kwargs)
    query = query.order_by(models.AMP.accession)
    data = query.offset(page * page_size).limit(page_size).all()
    for amp_obj in data:
        amp_obj.secondary_structure = None
        amp_obj.metadata = None
        amp_obj.num_genes = get_number_genes(amp_obj.accession, db)
    return mk_result(data, query.count(), page=page, page_size=page_size)


_number_genes_cache = {}
def get_number_genes(accession, db):
    global _number_genes_cache
    if r := _number_genes_cache.get(accession):
        return r
    r = db.query(func.count(models.GMSCMetadata.GMSC_accession)) \
                .where(models.GMSCMetadata.AMP == accession) \
                .scalar()
    _number_genes_cache[accession] = r
    return r


def get_amp(accession: str, db: Session):
    try:
        amp_obj = db.query(models.AMP).filter(models.AMP.accession == accession).one()
    except sqlalchemy.exc.NoResultFound:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    amp_obj.metadata = get_amp_metadata(amp_obj.accession, db, page=0, page_size=5)
    amp_obj.secondary_structure = utils.get_secondary_structure(amp_obj.sequence)
    return amp_obj


def get_amp_metadata(accession: str, db: Session, page: int, page_size: int):
    query = db.query(models.GMSCMetadata).filter(models.GMSCMetadata.AMP == accession)
    data = query.offset(page * page_size).limit(page_size).all()
    if len(data) == 0:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    return mk_result(data, total_items=query.count(), page_size=page_size, page=page)

def mk_result(data, total_items, page_size, page):
    return types.SimpleNamespace(
            data=data,
            info={
                'currentPage': page,
                'pageSize': page_size,
                'totalPage': math.ceil(total_items / page_size),
                'totalItem': total_items,
            })


def get_families(db: Session, page: int, page_size: int, **kwargs):
    query = db.query(distinct(models.AMP.family)).outerjoin(models.GMSCMetadata)
    # Mapping from filter keys to table columns
    metadata_cols = {
        'habitat': 'general_envo_name',
        'microbial_source': 'microbial_source',
        'sample': 'sample'}
    for key, value in kwargs.items():
        if value:
            query = query.filter(getattr(models.GMSCMetadata, metadata_cols[key]) == value)
    accessions = query.offset(page * page_size).limit(page_size).all()
    # if len(accessions) == 0:
    #     raise HTTPException(status_code=400, detail='invalid filter applied.')
    data = [get_family(accession, db=db) for accession, in accessions]
    return mk_result(data, query.count(), page=page, page_size=page_size)


def get_family(accession: str, db: Session):
    return dict(
        accession=accession,
        consensus_sequence=utils.cal_consensus_seq(accession),
        num_amps=db.query(func.count(models.AMP.accession).filter(models.AMP.family == accession)).scalar(),
        feature_statistics=get_fam_features(accession, db),
        distributions=get_distributions(accession, db),
        associated_amps=get_associated_amps(accession, db),
        downloads=get_fam_downloads(accession, db=db),
    )


def get_fam_metadata(accession: str, db: Session, page: int, page_size: int):
    amp_accessions = db.query(models.AMP.accession).filter(models.AMP.family == accession).all()
    amp_accessions = [accession for accession, in amp_accessions]
    # TODO FIX HERE
    m = db.query(models.GMSCMetadata).filter(models.GMSCMetadata.AMPSphere_code.in_(amp_accessions)). \
        offset(page * page_size).limit(page_size).all()
    return [row.__dict__ for row in m]

def get_amp_features(accession: str, db: Session):
    try:
        [seq] = db.query(models.AMP.sequence).filter(models.AMP.accession == accession).one()
    except sqlalchemy.exc.NoResultFound:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    return utils.get_amp_features(seq)

def get_fam_features(accession: str, db: Session):
    amps = db.query(models.AMP).filter(models.AMP.family == accession).all()
    features = [utils.get_amp_features(amp.sequence) for amp in amps]
    accessions = [amp.accession for amp in amps]
    if len(features) == 0:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    else:
        statistics = utils.recursive_round3(features)
        return {ix: row for ix, row in zip(accessions, statistics)}


def get_associated_amps(accession, db):
    amp_accessions = db.query(models.AMP.accession).filter(models.AMP.family == accession).all()
    return [accession for accession, in amp_accessions]


def get_distributions(accession: str, db: Session):
    if accession.startswith('AMP'):
        raw_data = db.query(models.GMSCMetadata).filter(models.GMSCMetadata.AMP == accession).all()
    elif accession.startswith('SPHERE'):
        raw_data = db.query(models.GMSCMetadata).outerjoin(models.AMP).filter(models.AMP.family == accession).all()
    else:
        raw_data = []
    if len(raw_data) == 0:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    return utils.compute_distribution_from_query_data(raw_data)


def get_fam_downloads(accession, db: Session):
    # TODO change prefix here for easier maintenance.
    q = db.query(models.AMP.family).filter(models.AMP.family == accession).first()
    in_db = bool(q)
    if not in_db:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    else:
        # FIXME Reimplement this.
        url_prefix = '{}/v1/families/{}/downloads/{}'.format(
            'https://ampsphere-api.big-data-biology.org',
            accession,
            accession)
        # FIX bug reported in issue 34, due to nginx proxy.
        print(url_prefix)
    path_bases = dict(
        alignment=str(url_prefix + '.aln'),
        sequences=str(url_prefix + '.faa'),
        # hmm_logo=str(url_prefix + '.png'),
        hmm_profile=str(url_prefix + '.hmm'),
        # sequence_logo=str(url_prefix + '.pdf'),
        tree_figure=str(url_prefix + '.ascii'),
        tree_nwk=str(url_prefix + '.nwk')
    )
    return {key: item.format(accession) for key, item in path_bases.items()}


def search_by_text(db: Session, text: str, page: int, page_size: int):
    """
    FIXME.
    :param query:
    :param db:
    :return:
    """
    # TODO retrieve text search result using the command: sqlite-utils search ampsphere_main_db/AMPSphere_v.2021-03.sqlite Metadata {search text}
    query = db.query(distinct(models.AMP.accession)).outerjoin(models.GMSCMetadata)
    # Consider blank space as + for query text.
    query = query.filter(or_(
        models.AMP.accession.like(text),
        models.AMP.family.like(text),
        models.GMSCMetadata.GMSC.like(text),
        models.GMSCMetadata.sample.like(text),
        models.GMSCMetadata.general_envo_name.like(text),
        models.GMSCMetadata.microbial_source.like(text),
    ))

    accessions = query.offset(page * page_size).limit(page_size).all()
    amps_data = [get_amp(accession, db) for accession, in accessions]
    return mk_result(amps_data, total_items=query.count(), page=page, page_size=page_size)


def get_statistics(db: Session):
    # TODO SPHEREs with num_amps < 8 should not be treated as families.
    # TODO display two numbers for num_genomes / num_metagenomes
    #               (analyzed_genomes..., num_...containing amps)
    # TODO FIX here.
    stats = db.query(models.Statistics).one()
    return dict(
        num_genes=stats.gmsc,
        num_amps=stats.amp,
        num_families=stats.family,
        num_habitats=stats.general_envo_name,
        num_genomes=db.query(func.count(distinct(models.GMSCMetadata.sample))).filter(
            models.GMSCMetadata.is_metagenomic == "False").scalar(),
        num_metagenomes=db.query(func.count(distinct(models.GMSCMetadata.sample))).filter(
            models.GMSCMetadata.is_metagenomic == "True").scalar(),
    )

def _all_used_taxa(db: Session):
    from collections import Counter
    used = set()
    for rank in 'pcofgs':
        used.update(db.execute(select(distinct(getattr(models.GMSCMetadata, f'microbial_source_{rank}')))).scalars().all())
    used.remove('')
    # If a genus is present, but all its elements are from the same species, remove it
    genus_species = db.query(models.GMSCMetadata.microbial_source_g, models.GMSCMetadata.microbial_source_s).distinct().all()
    cs = Counter(g for g, _ in genus_species)
    used -= set([k for k,v in cs.items() if v == 1])
    used = tuple(sorted(used))
    return used


_all_options = None
def get_all_options(db: Session):
    global _all_options
    if _all_options is None:
        habitat, = zip(*db.query(models.GMSCMetadata.general_envo_name).distinct())
        quality, = zip(*db.query(models.AMP.RNAcode).distinct())
        peplen_min, peplen_max, mw_min, mw_max, \
        pI_min, pI_max, charge_min, charge_max = db.query(
            func.min(models.AMP.length),
            func.max(models.AMP.length),
            func.min(models.AMP.molecular_weight),
            func.max(models.AMP.molecular_weight),
            func.min(models.AMP.isoelectric_point),
            func.max(models.AMP.isoelectric_point),
            func.min(models.AMP.charge),
            func.max(models.AMP.charge),
        ).first()
        round_floor = lambda x: Decimal(x).quantize(Decimal("0."), rounding=ROUND_FLOOR)
        round_ceiling = lambda x: Decimal(x).quantize(Decimal("0."), rounding=ROUND_CEILING)
        _all_options = dict(
            quality=quality,
            habitat=habitat,
            microbial_source=_all_used_taxa(db),
            pep_length=dict(min=int(peplen_min), max=int(peplen_max) + 1),
            molecular_weight=dict(min=round_floor(mw_min), max=round_ceiling(mw_max)),
            isoelectric_point=dict(min=round_floor(pI_min), max=round_ceiling(pI_max)),
            charge_at_pH_7=dict(min=round_floor(charge_min), max=round_ceiling(charge_max))
        )
    return _all_options


def entity_in_db(db, entity_type, accession):
    if entity_type == 'family':
        exists = db.query(db.query(models.AMP).filter(models.AMP.family == accession).exists()).scalar()
    elif entity_type in {'sample', 'genome'}:
        exists = db.query(db.query(models.GMSCMetadata).filter(models.GMSCMetadata.sample == accession).exists()).scalar()
    else:
        exists = False
    return exists


def mmseqs_search(seq: str, db):
    query_id = str(utils.uuid.uuid4())
    query_time_now = utils.datetime.now()
    tmp_dir = pathlib.Path(utils.cfg['tmp_dir'])
    input_seq_file = tmp_dir.joinpath(query_id + '.input')
    output_file = tmp_dir.joinpath(query_id + '.output')
    stdout_file = tmp_dir.joinpath(query_id + '.stdout')
    output_format = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qseq,tseq'
    # TODO calculate alignment string based on qaln,taln,gapopen,qstart,qend,tstart,alnlen
    # TODO add hint when the length of input sequence is not between 8 and 98
    # HINT: The search result may not reflect the reality as your sequence is too long/short.
    if not tmp_dir.exists():
        tmp_dir.mkdir(parents=True)
    with open(input_seq_file, 'w') as f:
        f.write(seq)
    # sensitivity = 1
    command_base = 'mmseqs createdb {query_seq} {query_seq}.mmseqsdb && ' \
                   'mmseqs search {query_seq}.mmseqsdb  {database} {out}.mmseqsdb {tmp_dir} -a && ' \
                   'mmseqs convertalis {query_seq}.mmseqsdb {database} {out}.mmseqsdb {out} --format-output {output_format}'
    command = command_base.format_map({
        'query_seq': input_seq_file,
        'database': utils.cfg['mmseqs_db'],
        'out': output_file,
        'tmp_dir': str(tmp_dir),
        'output_format': output_format,
        # 's': sensitivity
    })
    try:
        # TODO redirect the stdout to a temporary file and return its content when there is no match.
        with open(stdout_file, 'w') as f:
            utils.subprocess.run(command, shell=True, check=True, stdout=f)  ## FIXME
    except utils.subprocess.CalledProcessError as e:
        print('error when executing the command (code {})'.format(e))
        print(e.output)
        return None  # TODO better handle this
    else:
        columns = ['query_identifier', 'target_identifier', 'sequence_identity', 'alignment_length',
                   'number_mismatches', 'number_gap_openings', 'domain_start_position_query',
                   'domain_end_position_query', 'domain_start_position_target',
                   'domain_end_position_target', 'E_value', 'bit_score', 'seq_query', 'seq_target']
        try:
            df = pd.read_table(output_file, sep='\t', header=None)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=columns)
        df.columns = columns
        format_alignment0 = lambda x: format_alignment(
            x['seq_query'], x['seq_target'], x['bit_score'], x['domain_start_position_target'] - 1,
            x['domain_end_position_target']
        ).split('\n')[0:3]
        if df.shape[0] > 0:
            df['alignment_strings'] = df[['seq_query', 'seq_target', 'bit_score', 'domain_start_position_target',
                                          'domain_end_position_target']].apply(format_alignment0, axis=1)
            df['family'] = df['target_identifier'].apply(lambda x: get_family_by_amp(x, db))
        records = df.to_dict(orient='records')
        return records


def get_family_by_amp(amp_accession, db):
    family, = db.query(models.AMP.family).filter(models.AMP.accession == amp_accession).first()
    return family


def hmmscan_search(seq: str):
    # FIXME this doesn't work but reports no error, what happened?
    query_id = str(utils.uuid.uuid4())
    query_time_now = utils.datetime.now()
    tmp_dir = pathlib.Path(utils.cfg['tmp_dir'])
    input_seq_file = tmp_dir.joinpath(query_id + '.input')
    output_file = tmp_dir.joinpath(query_id + '.output')
    stdout_file = tmp_dir.joinpath(query_id + '.stdout')
    # TODO add hint when the length of input sequence is not between 8 and 98
    # HINT: The search result may not reflect the reality as your sequence is too long/short.

    if not tmp_dir.exists():
        tmp_dir.mkdir(parents=True)
    with open(input_seq_file, 'w') as f:
        f.write(seq)  # already in fasta format

    command_base = 'hmmscan --domtblout {out} {hmm_profiles} {query_seq}'
    command = command_base.format_map({
        'out': output_file,
        'hmm_profiles': utils.cfg['hmmprofile_db'],
        'query_seq': input_seq_file,
        'out_tmp': tmp_dir.joinpath(query_id + '.tmp')
    })
    try:
        # TODO redirect the stdout to a temporary file and return its content when there is no match.
        with open(stdout_file, 'w') as f:
            utils.subprocess.run(command, shell=True, check=True, stdout=f)  ## FIXME
    except utils.subprocess.CalledProcessError as e:
        print('error when executing the command (code {})'.format(e))
        print(e.output)
        return None
    else:
        columns = [
            'target_name', 'target_accession', 'target_length', 'query_name',
            'query_accession', 'query_length', 'E_value', 'score', 'bias',
            'domain_index', 'num_domain', 'c_Evalue', 'i_Evalue', 'score',
            'bias', 'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env',
            'to_env', 'acc', 'description_of_target']
        try:
            df = pd.read_table(output_file, header=2, skipfooter=10, sep=r'\s+', engine='python')
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=columns)
        df.columns = columns
        records = df.to_dict(orient='records')
        return records
