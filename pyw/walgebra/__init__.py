"""
W-algebra module - Quantum Drinfeld-Sokolov reduction and W-algebra constructions.
"""

from .quantum_ds_reduction import (
    AdmissibleWeightCandidate,
    CharacterEvaluation,
    LinearNonDegeneracyCondition,
    PrincipalAdmissibleLevelInput,
    QuantumDSReductionResult,
    ReducedModuleData,
    ReductionData,
    is_nondegenerate_affine_dynkin,
    quantum_ds_reduction_workflow,
)

__all__ = [
    "PrincipalAdmissibleLevelInput",
    "ReductionData",
    "LinearNonDegeneracyCondition",
    "is_nondegenerate_affine_dynkin",
    "AdmissibleWeightCandidate",
    "ReducedModuleData",
    "CharacterEvaluation",
    "QuantumDSReductionResult",
    "quantum_ds_reduction_workflow",
]
