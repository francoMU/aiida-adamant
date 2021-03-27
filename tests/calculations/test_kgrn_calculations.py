""" Tests for calculations

"""
import os
from pathlib import Path

from aiida.orm import StructureData, List

TEST_DIR = Path(__file__).parent / "data"
TEST_DIR = TEST_DIR.as_posix()


def test_creation_of_input_file():
    from aiida.plugins import DataFactory
    from aiida.orm import SinglefileData

    # Create an alloy
    alat = 4.
    cell = [[alat, 0., 0., ],
            [0., alat, 0., ],
            [0., 0., alat, ],
            ]
    structure = StructureData(cell=cell)

    structure.append_atom(position=(0., 0., 0.),
                          symbols=['Fe','Al'],
                          weights=[0.4,0.6])


    structure.append_atom(position=(0., 0., 0.),
                          symbols=['Ti','Al'],
                          weights=[0.3,0.7])


    print(structure)






def test_process(adamant_code):
    """Test running a calculation
    note this does not test that the expected outputs are created of output
    parsing"""
    from aiida.plugins import DataFactory, CalculationFactory
    from aiida.engine import run
    from aiida.orm import SinglefileData

    # Prepare input parameters
    parameters = {
        'niter': 200
    }

    KgrnParamsData = DataFactory('adamant.kgrn_data')

    params = KgrnParamsData(kgrn=parameters)

    test_file = SinglefileData(
        file=os.path.join(TEST_DIR, 'testfile.txt'))

    alat = 4.
    cell = [[alat, 0., 0., ],
            [0., alat, 0., ],
            [0., 0., alat, ],
            ]
    structure = StructureData(cell=cell)
    structure.append_atom(position=(0., 0., 0.), symbols='Fe')

    # set up calculation
    inputs = {
        'code': adamant_code,
        'kgrn': {
            'structure': structure,
            'params': params,
            'shape_function': test_file,
            'madelung_matrix': test_file,
            'atom_cfg': test_file,
            'transfer_matrix': List(list=[1, 2])
        },
        'metadata': {
            'options': {
                'max_wallclock_seconds': 30,
                'resources': {
                    'num_machines': 1,
                    'num_mpiprocs_per_machine': 1
                }
            },
        },
    }

    result = run(CalculationFactory('adamant.kgrn_calculation'), **inputs)

    print(result)

    # computed_diff = result['adamant'].get_content()
