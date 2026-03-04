"""Tests for quantum DS reduction MVP workflow."""

from fractions import Fraction

import pytest
from sage.all import CC, I, QQ

from pyw.utils.theta_functions import sl2_boundary_character
from pyw.fractional.level import FractionalLevel
from pyw.walgebra import (
    LinearNonDegeneracyCondition,
    PrincipalAdmissibleLevelInput,
    is_nondegenerate_affine_dynkin,
    quantum_ds_reduction_workflow,
)


def test_ds_workflow_builds_reduction_data_and_survivors():
    """Non-principal partition should still run end-to-end with explicit survivors."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=3,
        partition=[3, 1],
        level_input=PrincipalAdmissibleLevelInput(p=5, u=2),
    )

    assert result.reduction_data.orbit.partition == [3, 1]
    assert result.reduction_data.h_f_dimension >= 0
    assert result.reduction_data.delta_l_size >= 0
    assert len(result.surviving_candidates) > 0
    assert all((c.cohomology_h0_dim + c.cohomology_h1_dim) > 0 for c in result.surviving_candidates)


def test_ds_central_charge_principal_sl2_formula():
    """Principal sl2 central charge matches explicit closed formula."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=1,
        partition=[2],
        level_input=PrincipalAdmissibleLevelInput(p=3, u=2),
    )

    # k = -2 + 3/2 = -1/2 -> c = 1 - 6*(k+1)^2/(k+2) = 0
    c = result.all_modules[0].central_charge
    assert c == Fraction(0, 1)


def test_ds_reduced_conformal_weight_explicit_value():
    """Check an explicit reduced conformal weight from MVP policy."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=1,
        partition=[2],
        level_input=PrincipalAdmissibleLevelInput(p=4, u=3),
    )

    module = next(m for m in result.all_modules if m.candidate.finite_dynkin_labels == (1,))
    assert module.conformal_weight == Fraction(1, 16)


def test_ds_character_supported_sl2_boundary_path():
    """For supported sl2 boundary case, character is evaluated via theta utility."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=1,
        partition=[2],
        level_input=PrincipalAdmissibleLevelInput(p=2, u=3),
    )

    tau = I / 10
    z = QQ(1) / QQ(5)
    evaluation = result.evaluate_character(module_index=0, tau=tau, z=z, num_terms=50)

    assert evaluation.supported is True

    expected = sl2_boundary_character(u=3, j=0, tau=tau, z=z, num_terms=50)
    rel_err = abs(CC(evaluation.value) - CC(expected)) / abs(CC(expected))
    assert rel_err < 1e-12


def test_ds_character_unsupported_case_reports_reason():
    """Unsupported character paths should return a clear machine-testable reason."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=2,
        partition=[2, 1],
        level_input=PrincipalAdmissibleLevelInput(p=4, u=3),
    )

    evaluation = result.evaluate_character(module_index=0, tau=I / 10, z=Fraction(1, 5))
    assert evaluation.supported is False
    assert "sl2" in evaluation.reason.lower()


def test_ds_rejects_mismatched_fractional_level_input():
    """Workflow should reject FractionalLevel whose Cartan type differs from request."""
    level = FractionalLevel(["A", 1, 1], p=4, u=3)

    with pytest.raises(ValueError, match="does not match requested workflow input"):
        quantum_ds_reduction_workflow(
            cartan_type="A",
            rank=2,
            partition=[2, 1],
            level_input=level,
        )


def test_linear_nondegeneracy_condition_evaluator():
    """Generic linear filter accepts/rejects affine Dynkin labels as expected."""
    cond = LinearNonDegeneracyCondition(
        coeffs=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),
        rhs=Fraction(0, 1),
        relation="neq",
    )

    assert is_nondegenerate_affine_dynkin((Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)), [cond])
    assert not is_nondegenerate_affine_dynkin(
        (Fraction(0, 1), Fraction(-1, 1), Fraction(2, 1)), [cond]
    )


def test_workflow_applies_nondegeneracy_conditions():
    """Workflow should enforce custom affine-Dynkin non-degeneracy constraints."""
    cond = LinearNonDegeneracyCondition(
        coeffs=(Fraction(1, 1), Fraction(0, 1)),
        rhs=Fraction(0, 1),
        relation="neq",
    )

    base = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=1,
        partition=[2],
        level_input=PrincipalAdmissibleLevelInput(p=4, u=3),
    )
    constrained = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=1,
        partition=[2],
        level_input=PrincipalAdmissibleLevelInput(p=4, u=3),
        nondegeneracy_conditions=[cond],
    )

    # Phase-1 survivor rule uses N=0 BRST H^0; nondegeneracy can only reduce survivors.
    assert len(base.surviving_candidates) >= len(constrained.surviving_candidates)
    assert all(c.survives_reduction for c in constrained.surviving_candidates)


def test_workflow_exposes_truncated_brs_cohomology_dimensions():
    """Candidate data includes phase-1 N=0 BRST cohomology dimensions."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=1,
        partition=[2],
        level_input=PrincipalAdmissibleLevelInput(p=4, u=3),
    )

    by_labels = {m.candidate.finite_dynkin_labels: m.candidate for m in result.all_modules}

    # vacuum candidate has trivial N=0 cohomology and does not survive
    c0 = by_labels[(0,)]
    assert c0.reduction_grade == 0
    assert c0.cohomology_h0_dim == 0
    assert c0.cohomology_h1_dim == 0
    assert c0.survives_reduction is False

    # non-vacuum candidate with label 1 survives via nontrivial H^0
    c1 = by_labels[(1,)]
    assert c1.reduction_grade == 2
    assert c1.cohomology_h0_dim == 1
    assert c1.cohomology_h1_dim == 1
    assert c1.survives_reduction is True


def test_phase15_brs_nilpotency_check_passes_on_nontrivial_block():
    """Phase 1.5 should pass internal nilpotency check d1*d0=0."""
    result = quantum_ds_reduction_workflow(
        cartan_type="A",
        rank=3,
        partition=[3, 1],
        level_input=PrincipalAdmissibleLevelInput(p=5, u=2),
    )

    # If nilpotency failed, workflow would raise ValueError before returning.
    assert len(result.all_candidates) > 0
    assert all(c.cohomology_h0_dim >= 0 for c in result.all_candidates)
    assert all(c.cohomology_h1_dim >= 0 for c in result.all_candidates)
