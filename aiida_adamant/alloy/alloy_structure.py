import json
from typing import Union, List, Sequence, Optional

import numpy as np
from monty.json import MontyDecoder, MontyEncoder
from pymatgen.core import Lattice, Structure, PeriodicSite, FloatWithUnit, Unit
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from aiida_adamant.alloy.alloy_composition import AlloyComposition


class AlloyStructure(Structure):

    def __init__(self,
                 lattice: Union[List, np.ndarray, Lattice],
                 species: Sequence[AlloyComposition],
                 coords: Sequence[Sequence[float]],
                 charge: float = None,
                 validate_proximity: bool = False,
                 to_unit_cell: bool = False,
                 coords_are_cartesian: bool = False,
                 site_properties: dict = None,
                 properties: dict = None,
                 ):

        structure = Structure(lattice, species, coords, charge,
                              validate_proximity, to_unit_cell,
                              coords_are_cartesian)

        _site_properties = {'site_index': [],
                            'neq_site_index': []}

        symmetric_structure = SpacegroupAnalyzer(
            structure).get_symmetrized_structure()

        site_index = 0
        nonequivalent_site_index = 0
        for eq_sites in symmetric_structure.equivalent_sites:
            nonequivalent_site_index += 1
            for _ in eq_sites:
                site_index += 1
                _site_properties['site_index'].append(site_index)
                _site_properties['neq_site_index'].append(
                    nonequivalent_site_index)

        if site_properties is not None:
            _site_properties.update(site_properties)

        super().__init__(lattice,
                         species,
                         coords,
                         charge,
                         validate_proximity,
                         to_unit_cell,
                         coords_are_cartesian,
                         _site_properties)

        self._properties = properties if properties is not None else {}

    def __getattr__(self, attribute):

        props = object.__getattribute__(self, "_properties")
        if attribute in props:
            return props[attribute]

        raise AttributeError(f"{attribute} is not "
                             f"part of {self.__class__.__name__}")

    @property
    def properties(self):
        return self._properties

    @properties.setter
    def properties(self, value: dict):
        if value is not None:
            self._properties = value

    def as_dict(self, verbosity: int = 1, fmt=None, **kwargs) -> dict:

        lattice_dict = self._lattice.as_dict(verbosity=verbosity)
        del lattice_dict["@module"]
        del lattice_dict["@class"]

        alloy_dict = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "charge": self._charge,
            "lattice": lattice_dict,
            "sites": [],
        }

        for site in self:
            site_dict = site.as_dict(verbosity=verbosity)
            del site_dict["lattice"]
            del site_dict["@module"]
            del site_dict["@class"]
            del site_dict["label"]

            species_dict = site.species.as_dict()

            site_dict["species"] = species_dict

            alloy_dict["sites"].append(site_dict)

        if self._properties:
            alloy_dict["properties"] = self._properties

        return alloy_dict

    @classmethod
    def from_dict(cls, d, fmt=None):

        lattice = Lattice.from_dict(d["lattice"])
        charge = d.get("charge", None)

        sites = []

        for site_dict in d["sites"]:

            alloy_composition = \
                AlloyComposition.from_dict(site_dict["species"])

            props = site_dict.get("properties", None)
            if props is not None:
                for key in props.keys():
                    props[key] = json.loads(
                        json.dumps(props[key],
                                   cls=MontyEncoder), cls=MontyDecoder)

            lattice = lattice if lattice else Lattice.from_dict(d["lattice"])

            periodic_site = PeriodicSite(alloy_composition,
                                         coords=site_dict["abc"],
                                         lattice=lattice,
                                         properties=props
                                         )
            sites.append(periodic_site)

        obj = cls.from_sites(sites, charge=charge)

        obj.properties = d.get("properties", None)

        return obj

    @property
    def wigner_seitz_radius(self):
        wigner_seitz_radius = np.cbrt(self.volume * 3
                                      / (np.pi * self.num_sites * 4))

        _wigner_seitz_radius = FloatWithUnit(wigner_seitz_radius,
                                             unit=Unit('ang')).to('bohr')

        return float(_wigner_seitz_radius)

    @wigner_seitz_radius.setter
    def wigner_seitz_radius(self, value):
        volume = value ** 3 * 4 / 3 * np.pi * self.num_sites

        _volume = float(FloatWithUnit(volume, unit=Unit('bohr^3')).to('ang^3'))

        self.scale_lattice(_volume)

    def copy(self,
             site_properties: Optional[dict] = None,
             sanitize=False,
             properties: Optional[dict] = None):

        structure = super().copy(site_properties, sanitize)

        structure.properties.update(self._properties)

        if properties is not None:
            structure.properties.update(properties)

        site_props = self.site_properties
        if site_properties:
            site_props.update(site_properties)

        return structure

    def __eq__(self, other):
        if other is self:
            return True
        if other is None:
            return False
        if len(self) != len(other):
            return False
        if self.lattice != other.lattice:
            return False
        for site in self:
            if site not in other:
                return False
        if self.properties != other.properties:
            return False
        return True
