"""
This module implements a AlloyComposition class to represent a alloy
for the calculation with a Green's function based method
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from typing import Sequence, List, Dict, Any, Union

from monty.json import MSONable
from pymatgen.core import Composition
from pymatgen.core.periodic_table import get_el_sp, Species, Element

from aiida_adamant.utils.typing import MagneticParamsLike, ScreeningParamsLike

__author__ = "Franco Moitzi"
__version__ = "0.3"
__maintainer__ = "Franco Moitzi"
__email__ = "franco.moitzi@mcl.at"
__date__ = "Feb 13, 2021"


@dataclass
class MagneticParams(MSONable):
    """
    Class representing the magnetic model of an alloy component

    """
    is_paramagnetic: bool = True
    magnetic_model: str = "N"
    init_mag_mom: float = 1.9


@dataclass
class ScreeningParams(MSONable):
    """
    Class representing the electrostatics of an alloy component

    """
    alpha: float = 0.70
    beta: float = 1.2


class AlloyComposition(Composition):  # pylint: disable=too-many-ancestors
    """
    Class representing a AlloyComposition

    """
    def __init__(self,
                 elements: Sequence[Union[str, Element, Species]],
                 concentrations: Sequence[float],
                 magnetic_params: Sequence[MagneticParamsLike] = None,
                 screening_params: Sequence[ScreeningParamsLike] = None):
        """

        :param elements: list of elements
        :param concentrations: list of concentrations
        :param magnetic_params: list of magnetic model
        :param screening_params: list of electrostatic model
        """

        if len(elements) != len(concentrations):
            raise ValueError("Elements and concentrations doesn't have the "
                             "same size")

        self._magnetic_params = OrderedDict()
        for i, element in enumerate(elements):

            if magnetic_params is None:
                ele_params = MagneticParams()
            elif not isinstance(magnetic_params[i], MagneticParams):
                ele_params = MagneticParams(**magnetic_params[i])
            else:
                ele_params = magnetic_params[i]

            self._magnetic_params[get_el_sp(element)] = ele_params

        self._screening_params = OrderedDict()
        for i, element in enumerate(elements):

            if screening_params is None:
                ele_params = ScreeningParams()
            elif not isinstance(screening_params[i], ScreeningParams):
                ele_params = ScreeningParams(**screening_params[i])
            else:
                ele_params = screening_params[i]

            self._screening_params[get_el_sp(element)] = ele_params

        _elements = [str(get_el_sp(e)) for e in elements]

        self._is_paramagnetic = {
            get_el_sp(element): mag_param.is_paramagnetic
            for element, mag_param in self._magnetic_params.items()
        }

        super().__init__(strict=True, **dict(zip(_elements, concentrations)))

    @property
    def is_paramagnetic(self) -> Dict[Element, bool]:
        """

        :return: dictionary which describes with element is paramagnetic
        """
        return self._is_paramagnetic

    @property
    def concentrations(self) -> List[float]:
        """

        :return: list of concentrations
        """
        return list(self.values())

    @property
    def screening_params(self) -> Dict:
        """

        :return: a Dict of screening_params
        """
        return self._screening_params

    @property
    def magnetic_params(self) -> Dict:
        """

        :return: a Dict of magnetic_params
        """
        return self._magnetic_params

    def as_dict(self) -> Dict[str, Any]:
        """

        :return: serializable dict
        """

        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "elements": self.elements,
            "concentrations": self.concentrations,
            "magnetic_params": list(self.magnetic_params.values()),
            "screening_params": list(self.screening_params.values())
        }

    @classmethod
    def from_dict(cls, d) -> AlloyComposition:
        """

        :param d: dictionary representation of alloy
        :return: instance of the alloy composition
        """

        magnetic_params = [
            MagneticParams.from_dict(sd) for sd in d['magnetic_params']
        ]

        screening_params = [
            ScreeningParams.from_dict(sd) for sd in d['screening_params']
        ]

        elements = [Element.from_dict(e) for e in d['elements']]

        concentrations = [float(c) for c in d['concentrations']]

        return AlloyComposition(elements, concentrations, magnetic_params,
                                screening_params)

    def __eq__(self, other):
        """
        Compare two alloy composition

        :param other: alloy composition
        :return: bool if the two composition are the same
        """

        #  elements with amounts < Composition.amount_tolerance don't show up
        #  in the elmap, so checking len enables us to only check one
        #  compositions elements
        if len(self) != len(other):
            return False
        for ele, concentration in self.items():
            if abs(concentration - other[ele]) > Composition.amount_tolerance:
                return False
            if self.screening_params[ele] != other.screening_params[ele]:
                return False
            if self.magnetic_params[ele] != other.magnetic_params[ele]:
                return False

        return True
