"""
Test suite for pyw.core.weyl_group module
"""

import pytest

from sage.all import WeylGroup

from pyw.core.weyl_group import AffineWeylGroup


class TestAffineWeylGroup:
    """Test AffineWeylGroup class."""

    def test_init_a2_affine(self):
        """Test initialization with A2 affine type."""
        wg = AffineWeylGroup(["A", 2, 1])
        assert wg.cartan_type == ["A", 2, 1]
        assert not wg.extended

    def test_init_extended(self):
        """Test initialization with extended=True."""
        wg = AffineWeylGroup(["A", 2, 1], extended=True)
        assert wg.extended

    def test_init_finite_type(self):
        """Test initialization with finite type."""
        wg = AffineWeylGroup(["A", 2])
        assert wg.cartan_type == ["A", 2]

    def test_init_tuple_type(self):
        """Test initialization with tuple type."""
        wg = AffineWeylGroup(("A", 2, 1))
        assert wg.cartan_type == ("A", 2, 1)

    def test_weyl_group_property(self):
        """Test weyl_group property."""
        wg = AffineWeylGroup(["A", 2])
        sage_wg = wg.weyl_group
        # Should be a SageMath WeylGroup object
        assert hasattr(sage_wg, "simple_reflections")
        assert hasattr(sage_wg, "generators")

    def test_simple_reflections(self):
        """Test simple_reflections property."""
        wg = AffineWeylGroup(["A", 2])
        reflections = wg.simple_reflections
        # Should have at least 2 simple reflections for A2
        assert len(reflections) >= 2

    def test_dot_action(self):
        """Test dot_action method: w.λ = w(λ + ρ) - ρ."""
        wg = AffineWeylGroup(["A", 2])
        ws = wg.weyl_group.domain()

        # Get the Weyl vector (rho)
        rho = ws.rho() if hasattr(ws, "rho") else None

        # Create a simple weight
        fundamental = ws.fundamental_weights()
        weight = fundamental[1]

        # Get a simple reflection
        s1 = wg.weyl_group.simple_reflection(1)

        # Apply dot action
        result = wg.dot_action(s1, weight, rho)

        # Result should be a weight
        assert result is not None

    def test_translation_operator_affine(self):
        """Test translation_operator for affine type."""
        wg = AffineWeylGroup(["A", 2, 1])
        ws = wg.weyl_group.domain()

        # Get a simple root for translation
        roots = ws.simple_roots()
        alpha = list(roots.values())[0]

        # Create translation operator
        t = wg.translation_operator(alpha)

        # Translation should be a Weyl group element
        assert t is not None

    def test_translation_operator_extended(self):
        """Test translation_operator for extended affine Weyl group."""
        wg = AffineWeylGroup(["A", 2, 1], extended=True)

        # Extended Weyl group should support more translations
        ws = wg.weyl_group.domain()
        roots = ws.simple_roots()
        alpha = list(roots.values())[0]

        t = wg.translation_operator(alpha)
        assert t is not None

    def test_cartan_type_property(self):
        """Test cartan_type property."""
        wg = AffineWeylGroup(["A", 2, 1])
        assert wg.cartan_type == ["A", 2, 1]

    def test_extended_property(self):
        """Test extended property."""
        wg_standard = AffineWeylGroup(["A", 2, 1])
        wg_extended = AffineWeylGroup(["A", 2, 1], extended=True)

        assert not wg_standard.extended
        assert wg_extended.extended

    def test_str_representation(self):
        """Test string representation."""
        wg = AffineWeylGroup(["A", 2, 1])
        s = str(wg)
        assert "WeylGroup" in s or "A" in s
