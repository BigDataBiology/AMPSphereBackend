from src import main
from fastapi.testclient import TestClient
from src import crud, database, router

def test_basic_lookup():
    client = TestClient(main.app)
    r = client.get("/v1/amps/AMP10.000_001")
    data = r.json()
    assert data['sequence'] == 'FFGIGQQEMTLEEIGDKFGLTRERVRQIKEKAIRRLRQSNRSKLLKSYLG'

