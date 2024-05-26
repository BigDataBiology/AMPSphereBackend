import math
import pathlib
import subprocess
import types
import uuid
from datetime import datetime

import pandas as pd
from fastapi import HTTPException
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
from Bio.pairwise2 import format_alignment

from src import utils, database, schemas


def _page(items, page, page_size):
    return items[page * page_size: (page + 1) * page_size]

def get_filtered_amps(**criteria):
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
    amps = database.amps
    subquery = None
    criteria = {key: value for key, value in criteria.items() if value}

    for filter, value in criteria.items():
        if filter in {'habitat', 'sample_genome'}:
            if subquery is None:
                subquery = database.gmsc_metadata
            col = {
                'habitat': 'general_envo_name',
                'sample_genome': 'sample'
            }[filter]
            subquery = subquery.query(f'{col} == @value')
        elif filter == 'microbial_source':
            if value not in database.gtdb_taxon_to_rank:
                raise HTTPException(status_code=400, detail='wrong taxonomy name provided.')
            rank = database.gtdb_taxon_to_rank[value]
            if subquery is None:
                subquery = database.gmsc_metadata
            subquery = subquery.query(f'{rank} == @value')

        elif filter == 'exp_evidence' and value == 'Passed':
            amps = amps.query('metaproteomes == "Passed" or metatranscriptomes == "Passed"')
        elif filter == 'exp_evidence' and value == 'Failed':
            amps = amps.query('metaproteomes == "Failed" and metatranscriptomes == "Failed"')
        elif filter in {'family', 'antifam', 'RNAcode', 'coordinates'}:
            col = {
                'family': 'family',
                'antifam': 'Antifam',
                'RNAcode': 'RNAcode',
                'coordinates': 'coordinates'
                }[filter]
            amps = amps.query(f'{col} == @value')
        elif filter in {'pep_length_interval', 'mw_interval', 'pI_interval', 'charge_interval'}:
            min_v, max_v = value.split(',')
            min_v = float(min_v)
            max_v = float(max_v)
            col = {
                    'pep_length_interval': 'length',
                    'mw_interval': 'molecular_weight',
                    'pI_interval': 'isoelectric_point',
                    'charge_interval': 'charge'
                    }[filter]
            amps = amps.query(f'{min_v} <= {col} <= {max_v}')
    if subquery is not None:
        selected = set(subquery['AMP'])
        amps = amps[amps.index.map(selected.__contains__)]
    return amps


def get_amps(page: int = 0, page_size: int = 20, **kwargs):
    query = get_filtered_amps(**kwargs)
    data = _page(query, page, page_size)
    data = data.reset_index('accession').to_dict('records')
    data_as_obj = []
    for amp_obj in data:
        amp_obj = schemas.AMP(**amp_obj)
        amp_obj.secondary_structure = None
        amp_obj.metadata = None
        amp_obj.num_genes = database.number_genes_per_amp.get(amp_obj.accession, 0)
        data_as_obj.append(amp_obj)
    return mk_result(data_as_obj, len(query), page=page, page_size=page_size)


def get_amp(accession: str):
    try:
        amp = database.amps.loc[accession]
    except KeyError:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    amp_obj = schemas.AMP(accession=accession, **amp.to_dict())
    amp_obj.metadata = get_amp_metadata(accession, page=0, page_size=5)
    amp_obj.secondary_structure = utils.get_secondary_structure(amp_obj.sequence)
    return amp_obj


def get_amp_metadata(accession: str,  page: int, page_size: int):
    query = database.amp2gmsc.get(accession, [])

    if len(query) == 0:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    query = database.gmsc_metadata.loc[query]
    data = _page(query, page, page_size)
    data = data.reset_index().rename(columns={'accession':'GMSC_accession'}).to_dict('records')
    data_obj = []
    for it in data:
        if pd.isna(it['latitude']): it['latitude'] = None
        if pd.isna(it['longitude']): it['longitude'] = None
        data_obj.append(schemas.Metadata(**it))
    return mk_result(data_obj, len(query), page_size=page_size, page=page)

def mk_result(data, total_items, page_size, page):
    return types.SimpleNamespace(
            data=data,
            info={
                'currentPage': page,
                'pageSize': page_size,
                'totalPage': math.ceil(total_items / page_size),
                'totalItem': total_items,
            })


def get_families(page: int, page_size: int, habitat=None, microbial_source=None, sample=None):
    gmsc_metadata = database.gmsc_metadata
    query = ''

    if habitat is not None:
        query + 'general_envo_name == @habitat'
    if microbial_source is not None:
        try:
            rank = database.gtdb_taxon_to_rank[microbial_source]
        except KeyError:
            gmsc_metadata = gmsc_metadata.iloc[:0]
        else:
            if query: query += ' and '
            query += f'{rank} == @microbial_source'
    if sample is not None:
        if query: query += ' and '
        query += f'sample == @sample'
    if query:
        gmsc_metadata = gmsc_metadata.query(query)
    sel_amps = set(gmsc_metadata['AMP'])
    sel_families = sorted(set(database.amps.loc[sorted(sel_amps)]['family']))
    accessions = _page(sel_families, page, page_size)
    data = [get_family(accession) for accession in accessions]
    return mk_result(data, len(sel_families), page=page, page_size=page_size)


def get_family(accession: str):
    return dict(
        accession=accession,
        consensus_sequence=utils.cal_consensus_seq(accession),
        num_amps=database.amps.eval('family == @accession').sum(),
        feature_statistics=get_fam_features(accession),
        distributions=get_distributions(accession),
        associated_amps=get_associated_amps(accession),
        downloads=get_fam_downloads(accession),
    )


def get_amp_features(accession: str):
    try:
        seq = database.amps.loc[accession, 'sequence']
    except KeyError:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    return utils.get_amp_features(seq)

def get_fam_features(accession: str):
    amps = database.amps.query('family == @accession')
    features = amps.sequence.map(utils.get_amp_features).to_list()
    accessions = amps.index
    if len(features) == 0:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    else:
        statistics = utils.recursive_round3(features)
        return {ix: row for ix, row in zip(accessions, statistics)}


def get_associated_amps(accession):
    return database.amps.query('family == @accession').index.to_list()


def get_distributions(accession: str):
    accession
    if accession.startswith('AMP'):
        raw_data = database.gmsc_metadata.loc[database.amp2gmsc.get(accession, [])]
    elif accession.startswith('SPHERE'):
        sel_amps = database.amps.query('family == @accession').index
        raw_data = database.gmsc_metadata[database.gmsc_metadata['AMP'].map(set(sel_amps).__contains__)]
    else:
        raw_data = []
    if len(raw_data) == 0:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    return utils.compute_distribution_from_query_data(raw_data)


def get_fam_downloads(accession):
    # TODO change prefix here for easier maintenance.
    in_db = database.amps.eval('family == @accession').any()
    if not in_db:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    else:
        BASE_URL = 'https://ampsphere-api.big-data-biology.org/v1'
        url_prefix = f'{BASE_URL}/families/{accession}/downloads/{accession}'
    return dict(
        alignment=(url_prefix + '.aln'),
        sequences=(url_prefix + '.faa'),
        hmm_profile=(url_prefix + '.hmm'),
        tree_figure=(url_prefix + '.ascii'),
        tree_nwk=(url_prefix + '.nwk')
    )



def get_statistics():
    from collections import Counter
    # Using counter is ~3x faster than using pandas' value_counts
    return {
            'num_genes': len(database.gmsc_metadata),
            'num_amps':  len(database.amps),
            'num_families': sum(c >= 8 for c in Counter(database.amps['family']).values()),
            'num_habitats': len(set(database.gmsc_metadata['general_envo_name'])),
            'num_genomes':     len(set(database.gmsc_metadata.query('~is_metagenomic')['sample'])),
            'num_metagenomes': len(set(database.gmsc_metadata.query(' is_metagenomic')['sample'])),
    }

def _all_used_taxa():
    used = set()
    for rank in 'pcofgs':
        used.update(database.gmsc_metadata[f'microbial_source_{rank}'].dropna())
    # If a genus is present, but all its elements are from the same species, remove it
    nr_sp_g = database.gmsc_metadata\
                [['microbial_source_g', 'microbial_source_s']]\
                .drop_duplicates()\
                .groupby('microbial_source_g')\
                .size()
    used -= set(nr_sp_g[nr_sp_g <= 1].index)
    used = tuple(sorted(used))
    return used



_all_options = None
def get_all_options():
    global _all_options
    if _all_options is None:
        habitat = list(database.gmsc_metadata['general_envo_name'].unique())
        quality = list(database.amps['RNAcode'].unique())
        peplen_min = database.amps.sequence.str.len().min()
        peplen_max = database.amps.sequence.str.len().max()

        mw_min = database.amps.molecular_weight.min()
        mw_max = database.amps.molecular_weight.max()

        pI_min = database.amps.isoelectric_point.min()
        pI_max = database.amps.isoelectric_point.max()

        charge_min = database.amps.charge.min()
        charge_max = database.amps.charge.max()

        round_floor = lambda x: Decimal(x).quantize(Decimal("0."), rounding=ROUND_FLOOR)
        round_ceiling = lambda x: Decimal(x).quantize(Decimal("0."), rounding=ROUND_CEILING)
        _all_options = dict(
            quality=quality,
            habitat=habitat,
            microbial_source=_all_used_taxa(),
            pep_length=dict(min=int(peplen_min), max=int(peplen_max) + 1),
            molecular_weight=dict(min=round_floor(mw_min), max=round_ceiling(mw_max)),
            isoelectric_point=dict(min=round_floor(pI_min), max=round_ceiling(pI_max)),
            charge_at_pH_7=dict(min=round_floor(charge_min), max=round_ceiling(charge_max))
        )
    return _all_options


def entity_in_db(entity_type, accession):
    if entity_type == 'family':
        return accession in database.amps['family'].values
    if entity_type in {'sample', 'genome'}:
        return accession in database.gmsc_metadata['sample'].values
    return False


def mmseqs_search(seq: str):
    query_id = str(uuid.uuid4())
    query_time_now = datetime.now()
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
            subprocess.run(command, shell=True, check=True, stdout=f)  ## FIXME
    except subprocess.CalledProcessError as e:
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
            df['family'] = df['target_identifier'].map(database.amps['family'])
        records = df.to_dict(orient='records')
        return records


def hmmscan_search(seq: str):
    # FIXME this doesn't work but reports no error, what happened?
    query_id = str(uuid.uuid4())
    query_time_now = datetime.now()
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
            subprocess.run(command, shell=True, check=True, stdout=f)  ## FIXME
    except subprocess.CalledProcessError as e:
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

def exact_sequence_match(seq: str):
    if len(seq) < 8 or len(seq) > 98:
        return None
    seq = seq.replace(' ', '')
    amp = database.amps.query('sequence == @seq')
    if len(amp) == 0 and seq[0] == 'M':
        seq = seq[1:]
        amp = database.amps.query('sequence == @seq')
    if amp.shape[0] == 0:
        return None
    else:
        return amp.index[0]

