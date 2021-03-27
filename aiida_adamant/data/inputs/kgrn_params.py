"""
Datatype and methods for the KGRN params
"""
import copy
import json
import operator
from pathlib import Path
from typing import Mapping, Optional

from aiida.orm import Dict
from monty.io import zopen
from pymatgen.util.typing import PathLike

cwd = Path(__file__).parent
with open(cwd / "DEFAULT_PARAMS.json") as params:
    DEFAULT_PARAMS = json.loads(params.read())


class KgrnParamsData(Dict):
    """
    KgrnParamsData(kgrn=None)

    AiiDA compatible node representing a KGRN parameter data object.

    The parameters a validate before. The only things that is done here

    """

    def __init__(self, **kwargs):
        dictionary = kwargs.pop('kgrn', None)

        if dictionary is None:
            super().__init__(**kwargs)
        else:
            _dictionary = self._check_params(dictionary)
            super().__init__(dict=_dictionary)

    @staticmethod
    def _check_params(dictionary: Mapping) -> Mapping:
        """
        Check params of dictionary

        :param dictionary:
        :return:
        """

        _dictionary = dict(sorted(DEFAULT_PARAMS.items(),
                                  key=operator.itemgetter(0)))

        for param, value in dictionary.items():

            if param.lower() in DEFAULT_PARAMS:
                _dictionary[param.lower()] = value
            else:
                raise ValueError(f'Parameter {param.lower()} '
                                 f'is not part of the kgrn input')

        return _dictionary



