"""
Core module - SageMath wrappers for root systems, Weyl groups, and weight spaces.
"""

from .affine_lie_algebra import (
    AffineLieAlgebra,
    scalar_product,
    weyl_reflection,
)
from .affine_weight import (
    AffineWeight,
    AffineWeightKey,
    affine_weight,
    affine_weight_key,
    from_dynkin_labels,
)
from .bruhat import BruhatOrder, CosetRepresentative, ParabolicSubgroup
from .character import (
    BoundedKLOrbit,
    IntegrableModuleCharacter,
    KazhdanLusztigCharacter,
    KazhdanLusztigData,
    KLNumeratorTerm,
)
from .hybrid_affine_character import (
    BoundedAffinePositiveRoot,
    BoundedAffineVermaPartition,
    BoundedWeightTree,
    HybridCharacterStats,
    KazhdanLusztigFreudenthalCharacter,
    KLCandidate,
    affine_simple_root_coordinates,
    bounded_affine_positive_roots,
    bounded_kl_candidates_above,
    build_bounded_weight_tree,
    classically_dominant_affine_weight,
    compute_bounded_freudenthal_multiplicities,
    freudenthal_denominator,
    is_in_affine_positive_root_cone,
    kl_coefficient_at_weight,
)
from .kazhdan_lusztig import KazhdanLusztigPolynomials
from .root_system import AffineRootSystem
from .weight_space import FractionalWeightSpace
from .weyl_group import AffineWeylGroup

__all__ = [
    "AffineRootSystem",
    "AffineWeylGroup",
    "FractionalWeightSpace",
    "AffineLieAlgebra",
    "scalar_product",
    "weyl_reflection",
    # Di Francesco notation for affine weights
    "AffineWeight",
    "AffineWeightKey",
    "affine_weight_key",
    "affine_weight",
    "from_dynkin_labels",
    # Bruhat order and Kazhdan-Lusztig
    "BruhatOrder",
    "ParabolicSubgroup",
    "CosetRepresentative",
    "KazhdanLusztigPolynomials",
    "KazhdanLusztigData",
    "KLNumeratorTerm",
    # Character computation
    "IntegrableModuleCharacter",
    "KazhdanLusztigCharacter",
    "BoundedKLOrbit",
    "BoundedAffinePositiveRoot",
    "BoundedAffineVermaPartition",
    "BoundedWeightTree",
    "HybridCharacterStats",
    "KazhdanLusztigFreudenthalCharacter",
    "KLCandidate",
    "affine_simple_root_coordinates",
    "bounded_affine_positive_roots",
    "build_bounded_weight_tree",
    "classically_dominant_affine_weight",
    "compute_bounded_freudenthal_multiplicities",
    "freudenthal_denominator",
    "is_in_affine_positive_root_cone",
    "bounded_kl_candidates_above",
    "kl_coefficient_at_weight",
]
