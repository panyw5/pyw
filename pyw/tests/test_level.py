"""
Test suite for pyw.fractional.level module
"""

from fractions import Fraction

import pytest

from pyw.fractional.level import FractionalLevel


class TestFractionalLevel:
    """Test FractionalLevel class."""

    def test_init_valid_a2_affine(self):
        """Test initialization with valid A2 affine type."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        assert level.cartan_type == ["A", 2, 1]
        assert level.p == 4
        assert level.u == 3

    def test_level_property(self):
        """Test level property returns -h^∨ + p/u."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        # For A2, h^∨ = 3, so k = -3 + 4/3 = -5/3
        assert level.level == Fraction(-5, 3)

    def test_k_plus_h_vee_property(self):
        """Test k_plus_h_vee property returns p/u."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        assert level.k_plus_h_vee == Fraction(4, 3)

    def test_h_vee_property(self):
        """Test h_vee (dual Coxeter number) property."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        # For A2, h^∨ = 3
        assert level.h_vee == 3

    def test_lacety_a_type(self):
        """Test lacety for A-type (should be 1)."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        assert level.lacety == 1

    def test_lacety_d_type(self):
        """Test lacety for D-type (should be 1)."""
        level = FractionalLevel(["D", 4, 1], p=7, u=2)
        assert level.lacety == 1

    def test_lacety_e_type(self):
        """Test lacety for E-type (should be 1)."""
        level = FractionalLevel(["E", 6, 1], p=12, u=5)
        assert level.lacety == 1

    def test_lacety_b_type(self):
        """Test lacety for B-type (should be 2)."""
        level = FractionalLevel(["B", 3, 1], p=5, u=3)
        assert level.lacety == 2

    def test_lacety_c_type(self):
        """Test lacety for C-type (should be 2)."""
        level = FractionalLevel(["C", 2, 1], p=5, u=3)
        assert level.lacety == 2

    def test_lacety_f_type(self):
        """Test lacety for F-type (should be 2)."""
        level = FractionalLevel(["F", 4, 1], p=9, u=5)
        assert level.lacety == 2

    def test_lacety_g_type(self):
        """Test lacety for G-type (should be 3)."""
        level = FractionalLevel(["G", 2, 1], p=4, u=5)
        assert level.lacety == 3

    def test_validate_p_u_coprime(self):
        """Test validation: p and u must be coprime."""
        with pytest.raises(ValueError, match="p.*u.*coprime"):
            FractionalLevel(["A", 2, 1], p=4, u=2)

    def test_validate_u_lacety_coprime(self):
        """Test validation: u and lacety must be coprime."""
        # A2 has lacety=1, so any u works
        # B2 has lacety=2, so u must be odd
        with pytest.raises(ValueError, match="u.*lacety.*coprime"):
            FractionalLevel(["B", 2, 1], p=5, u=2)

    def test_validate_p_ge_h_vee(self):
        """Test validation: p >= h^∨."""
        with pytest.raises(ValueError, match="p.*>=.*h"):
            FractionalLevel(["A", 2, 1], p=2, u=3)

    def test_validate_p_positive(self):
        """Test validation: p must be positive."""
        with pytest.raises(ValueError, match="p.*positive"):
            FractionalLevel(["A", 2, 1], p=0, u=3)

    def test_validate_u_positive(self):
        """Test validation: u must be positive."""
        with pytest.raises(ValueError, match="u.*positive"):
            FractionalLevel(["A", 2, 1], p=4, u=0)

    def test_tuple_cartan_type(self):
        """Test with tuple Cartan type."""
        level = FractionalLevel(("A", 2, 1), p=4, u=3)
        assert level.level == Fraction(-5, 3)

    def test_finite_type(self):
        """Test with finite (non-affine) type."""
        level = FractionalLevel(["A", 2], p=4, u=3)
        # Finite types still work
        assert level.level == Fraction(-5, 3)  # -3 + 4/3

    def test_str_representation(self):
        """Test string representation."""
        level = FractionalLevel(["A", 2, 1], p=4, u=3)
        s = str(level)
        assert "A2~" in s or "A" in s
        assert "4" in s
        assert "3" in s

    @pytest.mark.parametrize("kind,rank,hdual", [("B", 3, 5), ("D", 4, 6), ("E", 6, 12), ("E", 7, 18), ("E", 8, 30), ("F", 4, 9)])
    def test_dual_coxeter_numbers(self, kind, rank, hdual):
        for ct in ([kind, rank], [kind, rank, 1]):
            level = FractionalLevel(ct, p=hdual, u=1)
            assert level.h_vee == hdual
            assert level.lookup_dual_coxeter_number() == hdual

    def test_g2_forbidden_denominator(self):
        with pytest.raises(ValueError, match="lacety.*coprime"):
            FractionalLevel(["G", 2, 1], p=4, u=3)

    def test_noninteger_parameters(self):
        with pytest.raises(ValueError, match="integers"):
            FractionalLevel(["A", 2, 1], p=4.5, u=3)
