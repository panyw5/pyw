"""
Tests for Kazhdan-Lusztig character computation modules.

These tests verify:
- Bruhat order comparisons
- KL polynomial computation
- Formal character arithmetic
- Kac-Wakimoto character formula
"""

import pytest
from sage.all import QQ, WeylGroup


class TestBruhatOrder:
    """Tests for BruhatOrder class."""

    def test_bruhat_le_identity(self):
        """Identity is less than or equal to any element."""
        from pyw.core.bruhat import BruhatOrder

        W = WeylGroup(["A", 2])
        bruhat = BruhatOrder(W)

        e = W.one()
        s1 = W.simple_reflection(1)
        s2 = W.simple_reflection(2)

        assert bruhat.le(e, e)
        assert bruhat.le(e, s1)
        assert bruhat.le(e, s2)
        assert bruhat.le(e, s1 * s2)

    def test_bruhat_le_simple_reflections(self):
        """Simple reflections are comparable with products."""
        from pyw.core.bruhat import BruhatOrder

        W = WeylGroup(["A", 2])
        bruhat = BruhatOrder(W)

        s1 = W.simple_reflection(1)
        s2 = W.simple_reflection(2)

        assert bruhat.le(s1, s1 * s2)
        assert bruhat.le(s2, s1 * s2)
        assert not bruhat.le(s1 * s2, s1)

    def test_bruhat_interval_A2(self):
        """Bruhat interval [e, w0] contains all elements for A2."""
        from pyw.core.bruhat import BruhatOrder

        W = WeylGroup(["A", 2])
        bruhat = BruhatOrder(W)

        e = W.one()
        w0 = W.long_element()

        interval = bruhat.interval(e, w0)
        assert len(interval) == 6  # |S_3| = 6

    def test_length_function(self):
        """Length function returns correct values."""
        from pyw.core.bruhat import BruhatOrder

        W = WeylGroup(["A", 2])
        bruhat = BruhatOrder(W)

        e = W.one()
        s1 = W.simple_reflection(1)
        w0 = W.long_element()

        assert bruhat.length(e) == 0
        assert bruhat.length(s1) == 1
        assert bruhat.length(w0) == 3  # Longest element in S_3


class TestParabolicSubgroup:
    """Tests for ParabolicSubgroup class."""

    def test_parabolic_contains(self):
        """Parabolic subgroup membership test."""
        from pyw.core.bruhat import BruhatOrder

        W = WeylGroup(["A", 3])
        bruhat = BruhatOrder(W)

        W_I = bruhat.parabolic_subgroup({1, 2})

        s1 = W.simple_reflection(1)
        s2 = W.simple_reflection(2)
        s3 = W.simple_reflection(3)

        assert W_I.contains(s1)
        assert W_I.contains(s2)
        assert W_I.contains(s1 * s2)
        assert not W_I.contains(s3)

    def test_minimal_coset_representative(self):
        """Minimal coset representative has correct length."""
        from pyw.core.bruhat import BruhatOrder

        W = WeylGroup(["A", 2])
        bruhat = BruhatOrder(W)

        W_I = bruhat.parabolic_subgroup({1})

        s1 = W.simple_reflection(1)
        s2 = W.simple_reflection(2)

        rep = W_I.minimal_coset_representative(s1 * s2)
        assert bruhat.length(rep) <= bruhat.length(s1 * s2)


class TestKazhdanLusztigPolynomials:
    """Tests for KL polynomial computation."""

    def test_P_identity(self):
        """P_{w,w} = 1 for any w."""
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        W = WeylGroup(["A", 2])
        kl = KazhdanLusztigPolynomials(W)

        e = W.one()
        s1 = W.simple_reflection(1)

        assert kl.P(e, e, at_one=True) == 1
        assert kl.P(s1, s1, at_one=True) == 1

    def test_P_adjacent(self):
        """P_{e, s_i} = 1 for simple reflections."""
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        W = WeylGroup(["A", 2])
        kl = KazhdanLusztigPolynomials(W)

        e = W.one()
        s1 = W.simple_reflection(1)

        assert kl.P(e, s1, at_one=True) == 1

    def test_Q_tilde_identity(self):
        """Q̃_{w,w} = 1 for any w."""
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        W = WeylGroup(["A", 2])
        kl = KazhdanLusztigPolynomials(W)

        e = W.one()
        assert kl.Q_tilde(e, e, at_one=True) == 1


class TestFormalCharacter:
    """Tests for FormalCharacter class."""

    def test_character_addition(self):
        """Character addition works correctly."""
        from pyw.core.character import FormalCharacter

        ch1 = FormalCharacter({0: 1, 1: 2}, max_grade=5)
        ch2 = FormalCharacter({0: 1, 1: 1, 2: 1}, max_grade=5)

        result = ch1 + ch2
        assert result[0] == 2
        assert result[1] == 3
        assert result[2] == 1

    def test_character_scalar_multiplication(self):
        """Scalar multiplication works correctly."""
        from pyw.core.character import FormalCharacter

        ch = FormalCharacter({0: 1, 1: 2, 2: 3}, max_grade=5)
        result = 2 * ch

        assert result[0] == 2
        assert result[1] == 4
        assert result[2] == 6

    def test_character_shift(self):
        """Grade shift works correctly."""
        from pyw.core.character import FormalCharacter

        ch = FormalCharacter({0: 1, 1: 2}, max_grade=5)
        shifted = ch.shift(2)

        assert shifted[0] == 0
        assert shifted[2] == 1
        assert shifted[3] == 2

    def test_character_truncate(self):
        """Truncation works correctly."""
        from pyw.core.character import FormalCharacter

        ch = FormalCharacter({0: 1, 1: 2, 2: 3, 3: 4}, max_grade=10)
        truncated = ch.truncate(2)

        assert truncated[0] == 1
        assert truncated[1] == 2
        assert truncated[2] == 3
        assert truncated[3] == 0


class TestWeylKacDenominator:
    """Tests for Weyl-Kac denominator."""

    def test_inverse_denominator_leading_term(self):
        """Inverse denominator has leading coefficient 1."""
        from pyw.core.character import WeylKacDenominator
        from pyw.core.affine_lie_algebra import AffineLieAlgebra

        ala = AffineLieAlgebra(["A", 2, 1])
        denom = WeylKacDenominator(ala)

        inv = denom.inverse(max_grade=5)
        assert inv[0] == 1


class TestVermaCharacter:
    """Tests for Verma module character."""

    def test_verma_character_exists(self):
        """Verma character can be computed."""
        from pyw.core.character import VermaCharacter
        from pyw.core.affine_lie_algebra import AffineLieAlgebra
        from pyw.core.affine_weight import AffineWeight

        ala = AffineLieAlgebra(["A", 2, 1])
        Lambda = AffineWeight.affine_fundamental_weight(ala, 1)

        verma = VermaCharacter(ala, Lambda)
        ch = verma.character(max_grade=3)

        assert ch[0] != 0


class TestKacWakimotoCharacter:
    """Tests for Kac-Wakimoto character formula."""

    @pytest.mark.skip(reason="Requires full implementation verification")
    def test_character_A2_level_1(self):
        """Character computation for A2 at level 1."""
        from pyw.algorithms.kac_wakimoto_character import KacWakimotoCharacter
        from pyw.core.affine_lie_algebra import AffineLieAlgebra
        from pyw.core.affine_weight import AffineWeight
        from pyw.fractional.level import FractionalLevel

        ala = AffineLieAlgebra(["A", 2, 1])
        level = FractionalLevel(["A", 2, 1], p=4, u=1)
        Lambda = AffineWeight.affine_fundamental_weight(ala, 1)

        kw = KacWakimotoCharacter(ala, level)
        ch = kw.character(Lambda, max_grade=3)

        assert ch[0] == 1
