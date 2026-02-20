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
from .bruhat import BruhatOrder, ParabolicSubgroup, CosetRepresentative
from .kazhdan_lusztig import KazhdanLusztigPolynomials
from .character import FormalCharacter, WeylKacDenominator, VermaCharacter

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
    # Bruhat order and Kazhdan-Lusztig
    "BruhatOrder",
    "ParabolicSubgroup",
    "CosetRepresentative",
    "KazhdanLusztigPolynomials",
    # Character computation
    "FormalCharacter",
    "WeylKacDenominator",
    "VermaCharacter",
]
