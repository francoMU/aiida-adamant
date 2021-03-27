from typing import Union, Dict, Any

from aiida_adamant.alloy.alloy_composition import MagneticParams, \
    ScreeningParams

MagneticParamsLike = Union[MagneticParams, Dict[str, Any]]
ScreeningParamsLike = Union[ScreeningParams, Dict[str, Any]]