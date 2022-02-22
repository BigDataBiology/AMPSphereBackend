import datatable as dt
import pandas as pd
import argparse
import tempfile
import pathlib
import subprocess


def main():
    args = get_args()
    check_dependencies()
    with tempfile.TemporaryDirectory(dir=args.tmp_dir) as tdir:
        extra_inputs_info = pd.read_table(args.extra_inputs_list, sep='\t')
        download_files(zenodo_id=args.zenodo_record_id, extra_inputs=extra_inputs_info, out=tdir)
        required_files = dict(   # should be
            aa_seq="AMPSphere_v.2021-03.faa.gz",
            family_align="AMPSphere_v.2021-03_families_alignment.tar.gz",
            family_tree="AMPSphere_v.2021-03_families_tree_nwk.tar.gz",
            na_seq="AMPSphere_v.2021-03.fna.xz",
            hosts="AMPSphere_v.2021-03.hosts.tsv.gz",
            locations="AMPSphere_v.2021-03.locations.tsv.gz",
            microontology="AMPSphere_v.2021-03.microontology.tsv.gz",
            origin_sample="AMPSphere_v.2021-03.origin_samples.tsv.gz",
            progenomes2_origins="AMPSphere_v.2021-03.species.tsv.gz",
            dramp_annotation="DRAMP_anno_AMPSphere_v.2021-03.parsed.tsv.gz",
            AMP_level_assessment="SPHERE_v.2021-03.levels_assessment.tsv.gz",
            quality_assessment='quality_assessment.tsv.xz',
            helical_wheels='AMPSphere_helicalwheels_latest_v2022.02.22.tgz'
        )
        check_downloaded_files(tdir, required_files)
        # TODO help needed.
        generate_main_database(required_files)
        generate_mmseqs_database(required_files)
        generate_hmmprofiles(required_files)
        generate_precomputed_data(required_files)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tmp_dir')
    parser.add_argument('-i', '--zenodo_id')
    parser.add_argument('--version-code')
    parser.add_argument('-e', '--extra_inputs_url', nargs='+')
    args = parser.parse_args()
    return args


def check_dependencies():
    return None


def download_files(zenodo_id, extra_inputs: pd.DataFrame, out):
    subprocess.Popen('zenodo_get -r {} -o {}'.format(zenodo_id, out))
    for file in extra_inputs.iterrows():
        filepath = out.join(file['filename'])
        MD5SUM = file['MD5SUM']
        subprocess.Popen('wget -c {url} -O | tee {filepath} | md5sum > {filepath}.MD5SUM'.format(url=file['url'], filepath=filepath))
        with open(filepath + '.MD5SUM', 'br') as f:
            assert f.read() == MD5SUM     # TODO add error information?
            # raise ValueError('MD5SUM checking failed, exiting...')


def check_downloaded_files(download_dir, expected_file_list):
    downloaded_files = pathlib.Path(download_dir).glob('*')
    for file, filename in expected_file_list.items():
        assert filename in downloaded_files  # TODO add error information?


