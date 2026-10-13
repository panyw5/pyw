"""
Test suite for pyw.fractional.admissible module
"""

from fractions import Fraction

import pytest
from sage.all import QQ, RootSystem

from pyw.core.weight_space import FractionalWeightSpace
from pyw.fractional.admissible import AdmissibleWeight
from pyw.fractional.level import FractionalLevel


def _reflection_coroots(cartan_type, depth):
    """Independent finite-depth Weyl orbit; no root-string formula is used."""
    lattice = RootSystem(cartan_type).coroot_lattice()
    reflections = lattice.simple_reflections()
    seen = set(lattice.simple_roots())
    frontier = seen.copy()
    for _ in range(depth):
        following = {s(b) for b in frontier for s in reflections}
        frontier = following - seen
        seen.update(frontier)
    return {b for b in seen if all(b.coefficient(i) >= 0 for i in lattice.index_set())}


class TestAdmissibleWeight:
    """Test AdmissibleWeight class."""

    def test_init_with_fractional_level(self):
        """Test initialization with FractionalLevel."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        assert aw.cartan_type == ("A", 2, 1)
        assert aw.level == level.level

    def test_init_with_numeric_level(self):
        """Test initialization with numeric level."""
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: 1, 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, 3)
        assert aw.level == Fraction(3)

    def test_setup_creates_weight_space(self):
        """Test _setup creates weight space."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        assert aw.weight_space is not None
        assert aw.root_system is not None

    def test_compute_affine_weyl_vector(self):
        """Test _compute_affine_weyl_vector."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        assert aw.rho_hat is not None

    def test_is_admissible_returns_bool(self):
        """Test is_admissible returns boolean."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        result = aw.is_admissible()
        assert isinstance(result, bool)

    def test_condition_1a_returns_bool(self):
        """Test _check_condition_1a returns boolean."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        result = aw._check_condition_1a()
        assert isinstance(result, bool)

    def test_condition_1b_returns_bool(self):
        """Test _check_condition_1b returns boolean."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        result = aw._check_condition_1b()
        assert isinstance(result, bool)

    def test_simple_admissible_weight(self):
        """Test with a simple admissible weight."""
        # For A2 affine with level k = -5/3, some weights are admissible
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        ws = FractionalWeightSpace(["A", 2, 1])

        # Try a dominant weight
        weight = ws.create_fractional_weight({0: Fraction(-5, 3)})

        aw = AdmissibleWeight(["A", 2, 1], weight, level)
        assert aw.is_admissible()

    def test_finite_type(self):
        """Test with finite (non-affine) type."""
        ws = FractionalWeightSpace(["A", 2])
        weight = ws.create_fractional_weight({1: 1, 2: 0})

        aw = AdmissibleWeight(["A", 2], weight, 1)
        assert aw.cartan_type == ("A", 2)

    @pytest.mark.parametrize("ct", [("A", 2), ("A", 2, 1)])
    def test_public_weyl_group_wrapper(self, ct):
        from pyw.core.weyl_group import AffineWeylGroup

        ws = FractionalWeightSpace(ct)
        aw = AdmissibleWeight(ct, ws.zero(), 0)
        wrapper = aw.weyl_group
        assert isinstance(wrapper, AffineWeylGroup)
        assert aw.weyl_group is wrapper
        assert tuple(wrapper.weyl_group.cartan_type()) == ct
        assert set(wrapper.simple_reflections) == set(aw.root_system.cartan_type().index_set())
        for reflection in wrapper.simple_reflections.values():
            assert reflection * reflection == wrapper.weyl_group.one()

    @pytest.mark.parametrize("ct", [("A", 2, 1), ("B", 3, 1), ("C", 3, 1), ("G", 2, 1)])
    def test_affine_coroot_strings_against_reflection_orbit(self, ct):
        # Non-vacuum rational labels make every finite component affect pairing.
        space = RootSystem(ct).weight_space(QQ, extended=True)
        indices = tuple(space.index_set())
        lambdas = space.fundamental_weights()
        weight = sum((QQ((2 * i + 1, i + 2)) * lambdas[i] for i in indices), space.zero())
        aw = AdmissibleWeight(ct, weight, weight.level())
        shifted = weight + space.rho()
        oracle = _reflection_coroots(ct, depth=32)
        oracle_vectors = {tuple(b.coefficient(i) for i in indices): b for b in oracle}
        strings = list(aw._affine_coroot_strings())
        assert len(strings) == len(RootSystem(ct[:2]).coroot_lattice().roots())

        generated = {}
        # Enough terms to cover every BFS root, independently bounded by its
        # affine-node coordinate. Every string has a positive affine step.
        max_affine = max(v[0] for v in oracle_vectors)
        for base, initial, step, direction in strings:
            assert direction[0] > 0
            for n in range(int(max_affine / direction[0]) + 1):
                vector = tuple(v + n * d for v, d in zip(base, direction, strict=True))
                assert vector not in generated
                generated[vector] = initial + n * step
                if n < 3:
                    # Both signs of each finite coroot, with the first three
                    # positive terms, must actually occur in the Weyl orbit.
                    assert vector in oracle_vectors, (ct, vector)
                    assert generated[vector] == shifted.scalar(oracle_vectors[vector])

        # Conversely, no positive real coroot found by reflections is missing.
        for vector, coroot in oracle_vectors.items():
            assert vector in generated, (ct, vector)
            assert generated[vector] == shifted.scalar(coroot)

    def test_tuple_cartan_type(self):
        """Test with tuple Cartan type."""
        level = FractionalLevel(("A", 2, 1), p=4, u=3)
        ws = FractionalWeightSpace(("A", 2, 1))
        weight = ws.create_fractional_weight({0: Fraction(-11, 3), 1: 1, 2: 1})

        aw = AdmissibleWeight(("A", 2, 1), weight, level)
        assert aw.cartan_type == ("A", 2, 1)

    def test_level_mismatch(self):
        ws = FractionalWeightSpace(["A", 2, 1])
        with pytest.raises(ValueError, match="central level"):
            AdmissibleWeight(["A", 2, 1], ws.zero(), Fraction(-5, 3))

    def test_finite_nonsimple_coroot_obstruction(self):
        ws = FractionalWeightSpace(["A", 2])
        weight = ws.create_fractional_weight({1: (-3, 2), 2: (-3, 2)})
        aw = AdmissibleWeight(["A", 2], weight, 0)
        # Simple shifted labels are -1/2, but their sum is -1.
        assert not aw._check_condition_1a()
        assert not aw.is_admissible()

    def test_finite_fractional_span_fails(self):
        ws = FractionalWeightSpace(["A", 2])
        aw = AdmissibleWeight(["A", 2], ws.create_fractional_weight({1: (1, 2), 2: (1, 2)}), 0)
        assert aw._check_condition_1a()
        assert not aw._check_condition_1b()

    @pytest.mark.parametrize("ct,p,u", [(["A", 2, 1], 4, 3), (["B", 3, 1], 5, 3), (["G", 2, 1], 4, 5)])
    def test_affine_vacuum_admissible(self, ct, p, u):
        level = FractionalLevel(ct, p, u)
        ws = FractionalWeightSpace(ct)
        aw = AdmissibleWeight(ct, ws.create_fractional_weight({0: level.level}), level)
        assert aw._check_condition_1a()
        assert aw._check_condition_1b()
        assert aw.is_admissible()

    def test_affine_nonsimple_coroot_obstruction(self):
        ws = FractionalWeightSpace(["A", 2, 1])
        weight = ws.create_fractional_weight({0: (7, 3), 1: (-3, 2), 2: (-3, 2)})
        aw = AdmissibleWeight(["A", 2, 1], weight, Fraction(-2, 3))
        assert not aw._check_condition_1a()

    def test_affine_critical_and_negative_shift(self):
        ws = FractionalWeightSpace(["A", 1, 1])
        for k in (-2, -3):
            aw = AdmissibleWeight(["A", 1, 1], ws.create_fractional_weight({0: k}), k)
            assert not aw.is_admissible()

    def test_twisted_type_explicitly_unsupported(self):
        ws = FractionalWeightSpace(["A", 2, 2])
        with pytest.raises(NotImplementedError, match="untwisted"):
            AdmissibleWeight(["A", 2, 2], ws.zero(), 0)

    def test_pyw_affine_weight_conversion(self):
        from sage.all import QQ

        from pyw.core.affine_lie_algebra import AffineLieAlgebra
        from pyw.core.affine_weight import from_dynkin_labels

        algebra = AffineLieAlgebra(["A", 2, 1])
        weight = from_dynkin_labels(algebra, {0: QQ((-5, 3)), 1: 0, 2: 0})
        assert AdmissibleWeight(["A", 2, 1], weight, Fraction(-5, 3)).is_admissible()
