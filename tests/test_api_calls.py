import json
import pandas as pd
from src import main
from fastapi.testclient import TestClient
from src import crud, database, router
import pytest

client = TestClient(main.app)

def test_root():
    response = client.get("/")
    assert response.status_code == 200

def test_basic_lookup():
    r = client.get("/v1/amps/AMP10.000_001")
    data = r.json()
    assert data['sequence'] == 'FFGIGQQEMTLEEIGDKFGLTRERVRQIKEKAIRRLRQSNRSKLLKSYLG'


def test_amps():
    r = client.get("v1/amps?page=0&page_size=5")
    assert json.load(open('tests/expected/amps5.json')) == r.json()


def test_sphere():
    r = client.get('/v1/families/SPHERE-III.000_264')
    assert r.json() == json.load(open('tests/expected/sphere.000_264.json', 'rt'))


def test_downloads():
    r = client.get("v1/downloads/GMSCMetadata.tsv")
    assert r.headers['content-type'].startswith('text/tab-separated-values')
    assert 'last-modified' in r.headers


def test_statistics():
    r = client.get("v1/statistics")
    assert r.json()['num_habitats'] == 73


def test_microbial_source_filter():
    data = client.get("v1/amps?page=0&page_size=100&microbial_source=Pseudoalteromonas%20luteoviolacea_F").json()
    assert len(data['data']) == 2
    assert [x['accession'] for x in data['data']] == ['AMP10.224_819', 'AMP10.723_664']



inputs = pd.read_table('tests/inputs.tsv', sep='\t').fillna('')
inputs['Test path'] = inputs.apply(
    lambda x: x['API path'].format(accession=x['Input accession']) +
              ('?{}'.format(x['Input query']) if x['Input query'] != '' else ''),
    axis=1)

@pytest.mark.parametrize('test_path, response_code', inputs[['Test path', 'Response code']].values.tolist())
def test_all(test_path, response_code):
    response = client.get(test_path)
    assert response.status_code == response_code
