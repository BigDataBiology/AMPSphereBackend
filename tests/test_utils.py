from src import utils
def test_get_secondary_structure():
    assert utils.get_secondary_structure("AAAAAAAA") == {'helix': 0.0, 'turn': 0.0, 'sheet': 1.0}
    assert utils.get_secondary_structure("AAAAAAAACCCCCCCCCCCC") == {'helix': 0.0, 'turn': 0.0, 'sheet': 0.4}
    assert utils.get_secondary_structure("AAAAAAAACCCCCCCCCCCCWWWWW") == {'helix': 0.2, 'turn': 0.0, 'sheet': 0.32}


def test_recursive_round3():
    assert utils.recursive_round3({'a': 1.222222555}) == {'a' : 1.222}
    assert utils.recursive_round3([{'a': 1.222222555}]) == [{'a' : 1.222}]
    assert utils.recursive_round3([{'a': {'b': 1.222222555}}]) == [{'a' : {'b': 1.222}}]
