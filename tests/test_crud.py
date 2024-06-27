from src import crud

def test_basic():
    accession = 'AMP10.123_123'
    amp = crud.get_amp(accession)
    assert amp.accession == accession

def test_filters():
    amps = crud.get_amps(page=0, page_size=10)
    assert len(amps.data) == 10

def test_get_all_options():
    options = crud.get_all_options()
    for k in ['pep_length', 'molecular_weight', 'charge_at_pH_7', 'isoelectric_point']:
        assert k in options
        assert options['pep_length']['min'] < options['pep_length']['max']
    assert len(options['microbial_source']) == 29504

def test_features():
    # This caused a failure at some point
    distributions = crud.get_distributions('SPHERE-III.000_007')
    assert 'lat' in distributions['geo'].keys()

