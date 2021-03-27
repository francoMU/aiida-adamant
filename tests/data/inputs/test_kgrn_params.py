def test_store_and_load_kgrn_data():
    from aiida.orm import load_node
    from aiida.plugins import DataFactory

    parameters = {
        'niter': 200
    }

    KgrnParamsData = DataFactory('adamant.kgrn_data')

    node = KgrnParamsData(kgrn=parameters)

    assert node.get_dict()['strt'] == 'A'
    assert not node.is_stored

    uuid = node.store().uuid

    assert node.is_stored

    node_loaded = load_node(uuid)

    assert node.get_dict() == node_loaded.get_dict()
