import json
from src import main
from fastapi.testclient import TestClient
from src import crud, database, router

client = TestClient(main.app)

def test_basic_lookup():
    r = client.get("/v1/amps/AMP10.000_001")
    data = r.json()
    assert data['sequence'] == 'FFGIGQQEMTLEEIGDKFGLTRERVRQIKEKAIRRLRQSNRSKLLKSYLG'

def test_amps():
    r = client.get("v1/amps?page=0&page_size=5")
    assert json.load(open('tests/expected/amps5.json')) == r.json()
