from collections import OrderedDict
from typing import Sequence, Union

import numpy as np

from monty.json import MSONable
from numpy.typing import ArrayLike
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.units import Energy, FloatWithUnit

from adamant_base.core.alloy.alloy_structure import AlloyStructure


class StructureEntries(MSONable):

    def __init__(self,
                 total_energy: Energy,
                 ):
        self._total_energy = total_energy

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "total_energy": str(self.total_energy),
        }

    @classmethod
    def from_dict(cls, d):
        return cls(FloatWithUnit.from_string(d['total_energy']))

    @property
    def total_energy(self) -> Energy:
        return self._total_energy

    def __str__(self):
        return str(self.total_energy)

    def __eq__(self, other):
        return self.total_energy == other.total_energy


class ComponentEntries(OrderedDict, MSONable):

    def __getitem__(self, component):
        component = get_el_sp(component)
        return super().__getitem__(component)


class ComponentEntry(MSONable):

    def __init__(self,
                 magnetic_moments: Union[float, ArrayLike]
                 ):

        self._magnetic_moments = np.atleast_1d(magnetic_moments)

        if len(self._magnetic_moments) == 2:
            self.is_paramagnetic = True
        else:
            self.is_paramagnetic = False

        self.total_magnetic_moment = float(np.sum(magnetic_moments))

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "magnetic_moments": self.magnetic_moments
        }

    @property
    def magnetic_moments(self) -> ArrayLike:
        return self._magnetic_moments

    @property
    def magnetic_splitting(self) -> float:

        if len(self.magnetic_moments) > 1:
            return float(np.diff(self._magnetic_moments))
        return 0.0

    def __str__(self):
        return str(self.magnetic_moments)

    def __eq__(self, other):
        return np.isclose(self.magnetic_moments, other.magnetic_moments)


class AlloyEntries(MSONable):

    def __init__(self,
                 alloy_structure: AlloyStructure,
                 structure_entries: StructureEntries,
                 site_entries: Sequence[Sequence[ComponentEntry]]
                 ):
        self._alloy_structure = alloy_structure
        self._structure_entries = structure_entries

        # ensure that everything is a real mutable list
        self._base_site_entries = [list(s) for s in site_entries]

        if len(site_entries) != alloy_structure.num_sites:
            raise ValueError("The length is not the same")

        self._site_entries = OrderedDict()
        for index_site, site in enumerate(self._alloy_structure):

            comp_mapping = ComponentEntries()

            site_entry = site_entries[index_site]

            for index_comp, comp in enumerate(site.species):
                comp_entry = site_entry[index_comp]

                comp_mapping[comp] = comp_entry

            self._site_entries[site] = comp_mapping

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "alloy_structure": self._alloy_structure.as_dict(),
            "structure_entries": self._structure_entries.as_dict(),
            "site_entries": self._base_site_entries,
        }

    @property
    def site_entries(self):
        return self._site_entries

    @property
    def total_energy(self):
        return self._structure_entries.total_energy

    @property
    def total_energy_per_site(self):
        return (self._structure_entries.total_energy /
                self._alloy_structure.num_sites)
