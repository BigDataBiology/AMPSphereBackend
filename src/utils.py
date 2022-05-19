import configparser
import pathlib
import subprocess

import numpy as np
from Bio.SeqUtils.ProtParam import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from datetime import datetime
from pprint import pprint
import uuid
import pandas as pd
from Bio import SearchIO, AlignIO
from Bio.Align import AlignInfo
from Bio.pairwise2 import format_alignment
import livingTree as lt
from collections import OrderedDict


def parse_config():
    parser = configparser.ConfigParser()
    parser.read('config/config.ini')
    return parser['DEFAULT']


cfg = parse_config()

# Protparam scales:
# kd → Kyte & Doolittle Index of Hydrophobicity
# Flex → Normalized average flexibility parameters (B-values)
# hw → Hopp & Wood Index of Hydrophilicity
# em → Emini Surface fractional probability (Surface Accessibility)

aalist = ['A', 'C', 'D', 'E',
          'F', 'G', 'H', 'I',
          'K', 'L', 'M', 'N',
          'P', 'Q', 'R', 'S',
          'T', 'V', 'Y', 'W']

# Colour scheme in Lesk (Introduction to Bioinformatics)
# Uses 5 groups (note Histidine): 
# Small nonpolar	    G, A, S, T                  Orange
# Hydrophobic	        C, V, I, L, P, F, Y, M, W   Green
# Polar                 N, Q, H		                Magenta
# Negatively charged    D, E		                Red
# Positively charged    K, R                        Blue
colorpallete = {'G': 'orange', 'A': 'orange', 'S': 'orange', 'T': 'orange',
                'C': 'g', 'V': 'g', 'I': 'g', 'L': 'g',
                'P': 'g', 'F': 'g', 'Y': 'g', 'M': 'g',
                'W': 'g', 'N': 'm', 'Q': 'm', 'H': 'm',
                'D': 'r', 'E': 'r', 'K': 'b', 'R': 'b'}

# these hydrophobicity scales are minmax organized
# the higher the value, more hydrophobic the aa is
scales = {'Parker': {'W': 1.0, 'F': 0.96, 'L': 0.96,
                     'M': 0.71, 'V': 0.685, 'Y': 0.595,
                     'C': 0.43, 'P': 0.395, 'A': 0.395,
                     'H': 0.395, 'R': 0.29, 'T': 0.24,
                     'G': 0.215, 'K': 0.215, 'Q': 0.2,
                     'S': 0.175, 'N': 0.15, 'E': 0.11,
                     'D': 0.0, 'I': 0.9},
          'ez': {'L': -4.92, 'I': -4.92, 'V': -4.04,
                 'F': -2.98, 'M': -2.35, 'W': -2.33,
                 'A': -1.81, 'C': -1.28, 'G': -0.94,
                 'Y': 0.14, 'T': 2.57, 'S': 3.4,
                 'H': 4.66, 'Q': 5.54, 'K': 5.55,
                 'N': 6.64, 'E': 6.81, 'D': 8.72,
                 'R': 14.92, 'P': 0.0}}


def get_amp_features(seq, include_graph_points=True):
    """
    :param seq:
    :return:
        Dict{
            MW: ...,
            Length: ...,
            Molar_extinction: [..., ...],
            Aromaticity: ...,
            GRAVY: ...,
            Instability_index: ...,
            Isoeletric_point: ...,
            Charget_at_pH_7: ...,
            Secondary_structure: [..., ..., ...],
        }
    """
    analyzed_seq = ProteinAnalysis(str(seq))

    out = {'Secondary_structure': dict(zip(['helix', 'turn', 'sheet'], analyzed_seq.secondary_structure_fraction())),
           'Length': analyzed_seq.length,
           'Molar_extinction': dict(zip(['cysteines_reduced', 'cystines_residues'],
                                        analyzed_seq.molar_extinction_coefficient())),
           'Aromaticity': round_3(analyzed_seq.aromaticity()),
           'GRAVY': round_3(analyzed_seq.gravy()),
           'MW': round_3(analyzed_seq.molecular_weight()),
           'Charge_at_pH_7': round_3(analyzed_seq.charge_at_pH(7.0)),
           'Instability_index': round_3(analyzed_seq.instability_index()),
           'Isoelectric_point': round_3(analyzed_seq.isoelectric_point()),
           'graph_points': get_graph_points(seq) if include_graph_points else None}
    return out


def get_secondary_structure(seq):
    analyzed_seq = ProteinAnalysis(str(seq))
    return dict(zip(['helix', 'turn', 'sheet'], analyzed_seq.secondary_structure_fraction()))


def get_graph_points(seq):
    graphs = ['transfer_energy', 'hydrophobicity_parker', 'surface_accessibility', 'flexibility']
    param_dict_list = [scales['ez'], scales['Parker'], ProtParamData.em, ProtParamData.Flex]
    wd = 5
    analyzed_seq = ProteinAnalysis(str(seq))
    graph_points = {}
    for graph, param_dict in zip(graphs, param_dict_list):
        val = analyzed_seq.protein_scale(window=wd, param_dict=param_dict)
        xval, xcolors = [], []
        for i in range(0, len(val)):
            xcolors.append(colorpallete[seq[i]])
            xval.append(f'{seq[i]}{i}')
        graph_points[graph] = {'type': 'line plot', 'x': xval, 'y': val, 'c': xcolors}
    return graph_points


def fam_download_file(accession: str, file: str):
    extention = file.split('.')[-1]
    # aln  fastas  hmm_logo  hmm_profiles  seqlogos  tree_figs  tree_nwk
    folders = dict(
        aln='aln',
        faa='fastas',
        png='hmm_logo',
        hmm='hmm_profiles',
        pdf='seqlogos',
        ascii='tree_figs',
        nwk='tree_nwk'
    )
    file = pathlib.Path(cfg['pre_computed_data']).joinpath('families').joinpath(folders[extention]).joinpath(file)
    # print('file exists? ', file.exists())
    return file


def get_downloads():
    return ['AMPSphere_latest.sqlite',
            'AMPSphere_latest.mmseqsdb',
            'AMPSphere_latest.hmm',
            'AMP.tsv',
            'GMSCMetadata.tsv']


def download(file_type: str):
    mapping = {
        'AMPSphere_latest.sqlite': './data/ampsphere_main_db/AMPSphere_latest.sqlite',
        'AMPSphere_latest.mmseqsdb': './data/mmseqs_db/AMPSphere_latest.mmseqsdb',
        'AMPSphere_latest.hmm': './data/hmmprofile_db/AMPSphere_latest.hmm',
        'AMP.tsv': 'data/tables/AMP.tsv',
        'GMSCMetadata.tsv': 'data/tables/GMSCMetadata.tsv',
    }
    return mapping[file_type]


def cal_consensus_seq(accession):
    file = pathlib.Path(cfg.get('pre_computed_data')).joinpath('families/aln/{}.aln'.format(accession))
    if file.exists():
        alignment = AlignIO.read(file, 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        return str(summary_align.dumb_consensus())
    else:
        return ''


def obj_to_dict(obj):
    return {key: value for key, value in obj.__dict__.items() if key != '_sa_instance'}


def round_3(num: float):
    return round(num, 3)


def df_to_formatted_json(df, sep="."):
    """
    The opposite of json_normalize
    """
    result = []
    for idx, row in df.iterrows():
        parsed_row = {}
        for col_label, v in row.items():
            keys = col_label.split(".")

            current = parsed_row
            for i, k in enumerate(keys):
                if i == len(keys) - 1:
                    current[k] = v
                else:
                    if k not in current.keys():
                        current[k] = {}
                    current = current[k]
        # save
        result.append(parsed_row)
    return result


def get_sunburst_data(paths_values, sep: str = None) -> dict:
    """
    :param paths_values:
    :param sep: each path is a list if sep = None.
    :return:
    """
    # print(paths_values)
    paths_values.columns = ['path', 'value']
    paths = paths_values['path']
    if sep:
        paths = paths.str.strip(sep).str.split(sep)
    else:
        sep = ';'

    # paths to set up the prefix tree
    # pprint('before prefix')
    # pprint(paths.tolist())
    paths = paths.apply(lambda x: ['Unknown' if i == '' else i for i in x])
    paths = paths.apply(lambda x: [sep.join(x[0:i]) for i in range(1, len(x) + 1)])
    # pprint('after prefix')
    # pprint(paths.tolist())
    identifiers = paths.apply(lambda x: x[-1])
    tree = lt.SuperTree()
    tree.create_node(identifier='')
    tree.from_paths(paths)
    tree.init_nodes_data(0)
    values = dict(zip(identifiers, paths_values['value'].tolist()))
    # print(values)
    tree.fill_with(values)
    tree.update_values()

    # print(vars(tree.get_node('a')))
    def rm_prefix(identifier):
        return identifier.split(sep)[-1]

    id_parent_value = list(zip(*[(rm_prefix(nid), rm_prefix(tree.parent(nid).identifier), tree.get_node(nid).data)
                                 for nid in tree.expand_tree(mode=tree.WIDTH) if nid != '']))
    return dict(zip(['labels', 'parents', 'values'], id_parent_value))


def compute_distribution_from_query_data(query_data):
    if len(query_data) > 0:
        metadata = pd.DataFrame([obj.__dict__ for obj in query_data]).drop(columns='_sa_instance_state')
        # print(metadata)
        color_map = {}
        metadata['latitude'] = metadata['latitude'].replace('', np.nan).astype(float).round(1)
        metadata['longitude'] = metadata['longitude'].replace('', np.nan).astype(float).round(1)
        metadata['habitat_type'] = pd.Categorical(metadata['general_envo_name'].apply(lambda x: x.split(':')[0]))
        # print(metadata[['habitat_type', 'microontology']])
        # metadata['color'] = metadata['habitat_type'].map(color_map)
        data = dict(
            geo=metadata[['AMP', 'latitude', 'longitude', 'habitat_type']].
                groupby(['latitude', 'longitude', 'habitat_type'], as_index=False, observed=True).size(),
            habitat=metadata[['AMP', 'general_envo_name', 'habitat_type']].
                groupby(['general_envo_name', 'habitat_type'], as_index=False, observed=True).size(),
            microbial_source=metadata[['AMP', 'microbial_source_s']].
                groupby('microbial_source_s', as_index=False).size()
        )
        names = {'latitude': 'lat', 'longitude': 'lon', 'AMP': 'size'}
        data['geo'].rename(columns=names, inplace=True)
        data['geo'] = data['geo'].to_dict(orient='list')
        # FIXME hierarchical structure generation.
        data['habitat'] = data['habitat'][data['habitat'].general_envo_name != '']
        if data['habitat'].shape[0] > 0:
            data['habitat'] = get_sunburst_data(data['habitat'][['general_envo_name', 'size']], sep=':')
        else:
            data['habitat'] = dict(zip(['labels', 'parents', 'values'], [[], [], []]))

        def simplify(data):
            data = data.sort_values(by='size', ascending=False).set_index('microbial_source_s')['size']
            data = OrderedDict(zip(data.index.tolist(), data.tolist()))
            print(data)
            others_count = 0
            result = OrderedDict()
            for key, value in data.items():
                if key == '':
                    others_count += value
                elif len(result) < 10:
                    result[key] = value
                else:
                    others_count += value
            result['others *'] = others_count
            return result

        data['microbial_source'] = simplify(data['microbial_source'])
        data['microbial_source'] = dict(
            labels=list(data['microbial_source'].keys()),
            values=list(data['microbial_source'].values())
        )
        print(data['microbial_source'])
        return data
    else:
        empty_sunburst = dict(
            labels=[''],
            parents=[''],
            values=[0.0]
        )
        empty_bubblemap = dict(
            lat=[0.0],
            lon=[0.0],
            size=[0.0]
        )
        return dict(
            host=empty_sunburst,
            habitat=empty_sunburst,
            geo=empty_bubblemap
        )
