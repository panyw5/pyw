"""
Core module - SageMath wrappers for root systems, Weyl groups, and weight spaces.
"""

from .root_system import AffineRootSystem
from .weyl_group import AffineWeylGroup
from .weight_space import FractionalWeightSpace
from .affine_lie_algebra import (
    AffineLieAlgebra,
    scalar_product,
    get_marks,
    get_comarks,
    weyl_reflection,
)
from .affine_weight import (
    AffineWeight,
    affine_weight,
    from_dynkin_labels,
)

__all__ = [
    "AffineRootSystem",
    "AffineWeylGroup",
    "FractionalWeightSpace",
    "AffineLieAlgebra",
    "scalar_product",
    "get_marks",
    "get_comarks",
    "weyl_reflection",
    # Di Francesco notation for affine weights
    "AffineWeight",
    "affine_weight",
    "from_dynkin_labels",
]
