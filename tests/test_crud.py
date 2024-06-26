from src import crud

def test_basic():
    amp = crud.get_amp('AMP10.123_123')
    assert amp.accession == 'AMP10.123_123'

def test_filters():
    amps = crud.get_amps(page=0, page_size=10)
    assert len(amps.data) == 10

def test_get_all_options():
    options = crud.get_all_options()
    for k in ['pep_length', 'molecular_weight', 'charge_at_pH_7', 'isoelectric_point']:
        assert k in options
        assert options['pep_length']['min'] < options['pep_length']['max']
