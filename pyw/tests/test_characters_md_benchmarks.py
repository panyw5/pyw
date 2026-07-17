"""Benchmark character coefficients against pyw/tests/characters.md.

These cases check that KazhdanLusztigCharacter and
KazhdanLusztigFreudenthalCharacter recover the known relative-grade
flavored prefixes (without the overall q^{-c/24} prefactor).
"""

import pytest


def _symbolic_equal(left, right) -> bool:
    from sage.all import SR

    return bool((SR(left) - SR(right)).expand() == 0)


def _grade_dimensions_from_series(series, order: int) -> list[int]:
    from sage.all import SR, var

    q = var("q")
    series = SR(series)
    free = [v for v in series.variables() if v != q]
    unflavored = series.subs({v: 1 for v in free}) if free else series
    return [int(unflavored.coefficient(q, grade)) for grade in range(order + 1)]


def _hybrid_flavored_series(character_by_depth, algebra):
    from sage.all import QQ, SR, var

    weight_lattice = algebra.affine_weight_lattice_sage()
    simple_roots = weight_lattice.simple_roots()
    b_vars = {index: var(f"b{index}") for index in range(1, algebra.rank + 1)}
    q = var("q")
    total = SR(0)
    for depth, weights in character_by_depth.items():
        term = SR(0)
        for weight, multiplicity in weights.items():
            if multiplicity == 0:
                continue
            monomial = SR(1)
            for index in range(1, algebra.rank + 1):
                exponent = algebra.scalar_product(
                    weight,
                    algebra.from_sagemath(simple_roots[index]),
                )
                monomial *= b_vars[index] ** QQ(exponent)
            term += multiplicity * monomial
        total += term * q**depth
    return total


def _assert_series_prefix(actual, expected_by_grade: dict[int, object]) -> None:
    from sage.all import var

    q = var("q")
    for grade, expected in expected_by_grade.items():
        assert _symbolic_equal(actual.coefficient(q, grade), expected), (
            f"grade {grade}: got {actual.coefficient(q, grade)}, "
            f"expected {expected}"
        )


def _d4_level_minus_two_vacuum_expected():
    """Relative-grade flavored prefix from characters.md for (D4)_{-2} vacuum."""
    from sage.all import SR, var

    b1, b2, b3, b4 = var("b1 b2 b3 b4")
    grade1 = (
        4
        + b1**2 / b2
        + b2 / b1**2
        + b2 / b3**2
        + b2 * (1 + 1 / b4**2)
        + (b1 * (1 + b2) * (b2 + b3**2) * (b2 + b4**2)) / (b2**2 * b3 * b4)
        + ((1 + b2) * (b2 + b3**2) * (b2 + b4**2)) / (b1 * b2 * b3 * b4)
        + (1 + b3**2 + b4**2) / b2
    )
    grade2 = (
        1
        / (b1**4 * b2**4 * b3**4 * b4**4)
        * (
            b1**8 * b2**2 * b3**4 * b4**4
            + b2**6 * b3**4 * b4**4
            + b1**7 * b2 * (1 + b2) * b3**3 * (b2 + b3**2) * b4**3 * (b2 + b4**2)
            + b1 * b2**4 * (1 + b2) * b3**3 * (b2 + b3**2) * b4**3 * (b2 + b4**2)
            + b1**5
            * b2
            * (1 + b2)
            * b3
            * (b2 + b3**2)
            * b4
            * (b2 + b4**2)
            * (
                3 * b2 * b3**2 * b4**2
                + b3**2 * b4**2 * (1 + b3**2 + b4**2)
                + b2**2 * (b4**2 + b3**2 * (1 + b4**2))
            )
            + b1**3
            * b2**2
            * (1 + b2)
            * b3
            * (b2 + b3**2)
            * b4
            * (b2 + b4**2)
            * (
                3 * b2 * b3**2 * b4**2
                + b3**2 * b4**2 * (1 + b3**2 + b4**2)
                + b2**2 * (b4**2 + b3**2 * (1 + b4**2))
            )
            + b1**6
            * b3**2
            * b4**2
            * (
                b2**6
                + b3**4 * b4**4
                + b2**5 * (1 + b3**2 + b4**2)
                + b2**4 * (1 + b3**2 + b4**2) ** 2
                + b2 * b3**2 * b4**2 * (b4**2 + b3**2 * (1 + b4**2))
                + b2**2 * (b4**2 + b3**2 * (1 + b4**2)) ** 2
                + b2**3
                * (
                    b4**2
                    + b4**4
                    + b3**4 * (1 + b4**2)
                    + b3**2 * (1 + 6 * b4**2 + b4**4)
                )
            )
            + b1**2
            * b2**2
            * b3**2
            * b4**2
            * (
                b2**6
                + b3**4 * b4**4
                + b2**5 * (1 + b3**2 + b4**2)
                + b2**4 * (1 + b3**2 + b4**2) ** 2
                + b2 * b3**2 * b4**2 * (b4**2 + b3**2 * (1 + b4**2))
                + b2**2 * (b4**2 + b3**2 * (1 + b4**2)) ** 2
                + b2**3
                * (
                    b4**2
                    + b4**4
                    + b3**4 * (1 + b4**2)
                    + b3**2 * (1 + 6 * b4**2 + b4**4)
                )
            )
            + b1**4
            * b2
            * (
                b2**6 * b3**2 * b4**2
                + b3**6 * b4**6
                + b2 * b3**4 * b4**4 * (1 + b3**2 + b4**2) ** 2
                + b2**5 * (b4**2 + b3**2 * (1 + b4**2)) ** 2
                + b2**4
                * b3**2
                * b4**2
                * (1 + b3**4 + 6 * b4**2 + b4**4 + 6 * b3**2 * (1 + b4**2))
                + b2**2
                * b3**2
                * b4**2
                * (
                    b4**4
                    + 6 * b3**2 * (b4**2 + b4**4)
                    + b3**4 * (1 + 6 * b4**2 + b4**4)
                )
                + b2**3
                * b3**2
                * b4**2
                * (
                    2 * b3**4 * (1 + b4**2)
                    + 2 * (b4**2 + b4**4)
                    + b3**2 * (2 + 17 * b4**2 + 2 * b4**4)
                )
            )
        )
    )
    return {
        0: SR(1),
        1: grade1,
        2: grade2,
    }


@pytest.mark.sage
def test_kl_character_sl2_level1_matches_characters_md():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 1, 1])
    b1 = var("b1")
    vacuum = AffineWeight.affine_fundamental_weight(algebra, 0)
    omega1 = AffineWeight.affine_fundamental_weight(algebra, 1)
    engine = KazhdanLusztigCharacter(algebra)

    vacuum_expected = {
        0: SR(1),
        1: 1 + 1 / b1**2 + b1**2,
        2: 2 + 1 / b1**2 + b1**2,
        3: 3 + 2 / b1**2 + 2 * b1**2,
        4: 5 + 1 / b1**4 + 3 / b1**2 + 3 * b1**2 + b1**4,
    }
    omega1_expected = {
        0: 1 / b1 + b1,
        1: 1 / b1 + b1,
        2: 1 / b1**3 + 2 / b1 + 2 * b1 + b1**3,
        3: (1 + b1**2) ** 3 / b1**3,
    }

    vacuum_series = engine.character(vacuum, order=4)
    omega1_series = engine.character(omega1, order=3)

    _assert_series_prefix(vacuum_series, vacuum_expected)
    _assert_series_prefix(omega1_series, omega1_expected)
    assert _grade_dimensions_from_series(vacuum_series, 4) == [1, 3, 4, 7, 13]
    assert _grade_dimensions_from_series(omega1_series, 3) == [2, 2, 6, 8]


@pytest.mark.sage
def test_kl_character_sl2_level2_matches_characters_md():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 1, 1])
    b1 = var("b1")
    vacuum = 2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    omega1 = (
        AffineWeight.affine_fundamental_weight(algebra, 0)
        + AffineWeight.affine_fundamental_weight(algebra, 1)
    )
    two_omega1 = 2 * AffineWeight.affine_fundamental_weight(algebra, 1)
    engine = KazhdanLusztigCharacter(algebra)

    vacuum_expected = {
        0: SR(1),
        1: 1 + 1 / b1**2 + b1**2,
        2: (1 + b1**2 + b1**4) ** 2 / b1**4,
        3: 5 + 1 / b1**4 + 4 / b1**2 + 4 * b1**2 + b1**4,
    }
    omega1_expected = {
        0: 1 / b1 + b1,
        1: 1 / b1**3 + 2 / b1 + 2 * b1 + b1**3,
        2: 2 * (1 + 2 * b1**2 + 2 * b1**4 + b1**6) / b1**3,
        3: 1 / b1**5 + 4 / b1**3 + 8 / b1 + 8 * b1 + 4 * b1**3 + b1**5,
    }
    two_omega1_expected = {
        0: 1 + 1 / b1**2 + b1**2,
        1: 2 + 1 / b1**2 + b1**2,
        2: (1 + b1**2) ** 2 * (1 + b1**2 + b1**4) / b1**4,
        3: 7 + 2 / b1**4 + 5 / b1**2 + 5 * b1**2 + 2 * b1**4,
    }

    vacuum_series = engine.character(vacuum, order=3)
    omega1_series = engine.character(omega1, order=3)
    two_omega1_series = engine.character(two_omega1, order=3)

    _assert_series_prefix(vacuum_series, vacuum_expected)
    _assert_series_prefix(omega1_series, omega1_expected)
    _assert_series_prefix(two_omega1_series, two_omega1_expected)
    assert _grade_dimensions_from_series(vacuum_series, 3) == [1, 3, 9, 15]
    assert _grade_dimensions_from_series(omega1_series, 3) == [2, 6, 12, 26]
    assert _grade_dimensions_from_series(two_omega1_series, 3) == [3, 4, 12, 21]


@pytest.mark.sage
def test_kl_character_sl3_level1_matches_characters_md():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 2, 1])
    b1 = var("b1")
    b2 = var("b2")
    vacuum = AffineWeight.affine_fundamental_weight(algebra, 0)
    omega1 = AffineWeight.affine_fundamental_weight(algebra, 1)
    omega2 = AffineWeight.affine_fundamental_weight(algebra, 2)
    engine = KazhdanLusztigCharacter(algebra)

    vacuum_expected = {
        0: SR(1),
        1: (
            2
            + b1**2 / b2
            + b2 / b1**2
            + b1 * (1 / b2**2 + b2)
            + (1 / b2 + b2**2) / b1
        ),
        2: (
            5
            + (2 * b1**2) / b2
            + (2 * b2) / b1**2
            + (2 * b1 * (1 + b2**3)) / b2**2
            + (2 * (1 + b2**3)) / (b1 * b2)
        ),
    }
    omega1_expected = {
        0: b1 + 1 / b2 + b2 / b1,
        1: (b1**2 + b2 + b1 * b2**2) ** 2 / (b1**2 * b2**2),
    }
    omega2_expected = {
        0: 1 / b1 + b1 / b2 + b2,
        1: (b1 + b1**2 * b2 + b2**2) ** 2 / (b1**2 * b2**2),
    }

    vacuum_series = engine.character(vacuum, order=2)
    omega1_series = engine.character(omega1, order=1)
    omega2_series = engine.character(omega2, order=1)

    _assert_series_prefix(vacuum_series, vacuum_expected)
    _assert_series_prefix(omega1_series, omega1_expected)
    _assert_series_prefix(omega2_series, omega2_expected)
    assert _grade_dimensions_from_series(vacuum_series, 2) == [1, 8, 17]
    assert _grade_dimensions_from_series(omega1_series, 1) == [3, 9]
    assert _grade_dimensions_from_series(omega2_series, 1) == [3, 9]


@pytest.mark.sage
def test_hybrid_character_matches_characters_md_sl2_level1():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter

    algebra = AffineLieAlgebra(["A", 1, 1])
    b1 = var("b1")
    vacuum = AffineWeight.affine_fundamental_weight(algebra, 0)
    omega1 = AffineWeight.affine_fundamental_weight(algebra, 1)
    engine = KazhdanLusztigFreudenthalCharacter(algebra)

    vacuum_character = engine.character(vacuum, order=4)
    omega1_character = engine.character(omega1, order=3)
    vacuum_series = _hybrid_flavored_series(vacuum_character, algebra)
    omega1_series = _hybrid_flavored_series(omega1_character, algebra)

    vacuum_expected = {
        0: SR(1),
        1: 1 + 1 / b1**2 + b1**2,
        2: 2 + 1 / b1**2 + b1**2,
        3: 3 + 2 / b1**2 + 2 * b1**2,
        4: 5 + 1 / b1**4 + 3 / b1**2 + 3 * b1**2 + b1**4,
    }
    omega1_expected = {
        0: 1 / b1 + b1,
        1: 1 / b1 + b1,
        2: 1 / b1**3 + 2 / b1 + 2 * b1 + b1**3,
        3: (1 + b1**2) ** 3 / b1**3,
    }

    _assert_series_prefix(vacuum_series, vacuum_expected)
    _assert_series_prefix(omega1_series, omega1_expected)
    assert [sum(vacuum_character[depth].values()) for depth in range(5)] == [
        1,
        3,
        4,
        7,
        13,
    ]
    assert [sum(omega1_character[depth].values()) for depth in range(4)] == [
        2,
        2,
        6,
        8,
    ]


@pytest.mark.sage
def test_hybrid_character_matches_characters_md_sl2_level2_and_sl3():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter

    a1 = AffineLieAlgebra(["A", 1, 1])
    a2 = AffineLieAlgebra(["A", 2, 1])
    b1 = var("b1")
    b2 = var("b2")

    vacuum_a1 = 2 * AffineWeight.affine_fundamental_weight(a1, 0)
    vacuum_a2 = AffineWeight.affine_fundamental_weight(a2, 0)
    engine_a1 = KazhdanLusztigFreudenthalCharacter(a1)
    engine_a2 = KazhdanLusztigFreudenthalCharacter(a2)

    vacuum_a1_character = engine_a1.character(vacuum_a1, order=3)
    vacuum_a2_character = engine_a2.character(vacuum_a2, order=2)
    vacuum_a1_series = _hybrid_flavored_series(vacuum_a1_character, a1)
    vacuum_a2_series = _hybrid_flavored_series(vacuum_a2_character, a2)

    vacuum_a1_expected = {
        0: SR(1),
        1: 1 + 1 / b1**2 + b1**2,
        2: (1 + b1**2 + b1**4) ** 2 / b1**4,
        3: 5 + 1 / b1**4 + 4 / b1**2 + 4 * b1**2 + b1**4,
    }
    vacuum_a2_expected = {
        0: SR(1),
        1: (
            2
            + b1**2 / b2
            + b2 / b1**2
            + b1 * (1 / b2**2 + b2)
            + (1 / b2 + b2**2) / b1
        ),
        2: (
            5
            + (2 * b1**2) / b2
            + (2 * b2) / b1**2
            + (2 * b1 * (1 + b2**3)) / b2**2
            + (2 * (1 + b2**3)) / (b1 * b2)
        ),
    }

    _assert_series_prefix(vacuum_a1_series, vacuum_a1_expected)
    _assert_series_prefix(vacuum_a2_series, vacuum_a2_expected)
    assert [sum(vacuum_a1_character[depth].values()) for depth in range(4)] == [
        1,
        3,
        9,
        15,
    ]
    assert [sum(vacuum_a2_character[depth].values()) for depth in range(3)] == [
        1,
        8,
        17,
    ]


@pytest.mark.sage
def test_hybrid_and_kl_agree_on_characters_md_sl2_vacuum():
    from sage.all import var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter

    algebra = AffineLieAlgebra(["A", 1, 1])
    vacuum = AffineWeight.affine_fundamental_weight(algebra, 0)
    order = 3
    q = var("q")

    kl_series = KazhdanLusztigCharacter(algebra).character(vacuum, order=order)
    hybrid_series = _hybrid_flavored_series(
        KazhdanLusztigFreudenthalCharacter(algebra).character(vacuum, order=order),
        algebra,
    )

    for grade in range(order + 1):
        assert _symbolic_equal(
            kl_series.coefficient(q, grade),
            hybrid_series.coefficient(q, grade),
        )


@pytest.mark.sage
def test_hybrid_character_d4_level_minus_two_matches_characters_md():
    """(D4)_{-2} vacuum flavored prefix through order 2 from characters.md."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter

    algebra = AffineLieAlgebra(["D", 4, 1])
    vacuum = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    engine = KazhdanLusztigFreudenthalCharacter(algebra)

    character = engine.character(vacuum, order=2)
    series = _hybrid_flavored_series(character, algebra)
    expected = _d4_level_minus_two_vacuum_expected()

    _assert_series_prefix(series, expected)
    assert [sum(character[depth].values()) for depth in range(3)] == [1, 28, 329]


@pytest.mark.sage
def test_kl_character_d4_level_minus_two_matches_characters_md(tmp_path):
    """Eager KL recovers the characters.md (D4)_{-2} vacuum prefix through order 1."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

    algebra = AffineLieAlgebra(["D", 4, 1])
    vacuum = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    engine = KazhdanLusztigCharacter(algebra)
    engine.kl_polynomial = KazhdanLusztigPolynomials(
        algebra.affine_weyl_group_sage(),
        cache_dir=tmp_path,
        persistent_cache=False,
    )

    series = engine.character(vacuum, order=1)
    expected = _d4_level_minus_two_vacuum_expected()
    _assert_series_prefix(series, {0: expected[0], 1: expected[1]})
    assert _grade_dimensions_from_series(series, 1) == [1, 28]


@pytest.mark.sage
def test_kl_and_hybrid_reject_characters_md_boundary_weights():
    """Boundary cases in characters.md need the theta path, not KL/hybrid.

    The vacuum modules at k=-4/3 (sl2) and k=-3/2 (sl3) have fractional affine
    Dynkin labels after the dominant-Lambda search. Sage's ordinary affine weight
    lattice is integral, so the KL orbit preparation fails before any character
    coefficient is produced. Numerical characters.md checks for these modules live
    in test_theta_functions.py.
    """
    from sage.all import QQ

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter

    a1 = AffineLieAlgebra(["A", 1, 1])
    a2 = AffineLieAlgebra(["A", 2, 1])
    cases = [
        (a1, QQ(-4) / 3 * AffineWeight.affine_fundamental_weight(a1, 0)),
        (a2, QQ(-3) / 2 * AffineWeight.affine_fundamental_weight(a2, 0)),
    ]

    for algebra, highest_weight in cases:
        with pytest.raises(TypeError, match="rational to integer"):
            KazhdanLusztigCharacter(algebra).character(highest_weight, order=0)
        with pytest.raises(TypeError, match="rational to integer"):
            KazhdanLusztigFreudenthalCharacter(algebra).character(
                highest_weight, order=0
            )


@pytest.mark.sage
def test_characters_md_boundary_series_via_theta_path():
    """Cross-check characters.md boundary prefixes through the theta formulas."""
    from sage.all import CC, I

    from pyw.utils.theta_functions import (
        sl2_boundary_character,
        sl2_boundary_vacuum_character,
        sl3_boundary_vacuum_character,
    )

    from pyw.tests.test_theta_functions import (
        relative_error,
        relative_error_with_reference,
        sl2_boundary_nonvacuum_j1_leading_term,
        sl2_boundary_nonvacuum_j1_series_first_4_terms,
        sl2_boundary_nonvacuum_j2_leading_term,
        sl2_boundary_nonvacuum_j2_series_first_4_terms,
        sl2_vacuum_series_first_4_terms,
        sl3_vacuum_series_first_2_terms,
    )

    tau = 1.7 * I
    z = 0.23 * I
    z1 = 0.19 + 0.07 * I
    z2 = 0.13 + 0.04 * I

    vac = sl2_boundary_vacuum_character(3, tau, z, num_terms=500)
    assert relative_error(vac, sl2_vacuum_series_first_4_terms(tau, z)) < 1e-6

    j1 = sl2_boundary_character(3, 1, tau, z, num_terms=500)
    j1_series = sl2_boundary_nonvacuum_j1_series_first_4_terms(tau, z)
    j1_scale = CC(j1) / CC(sl2_boundary_nonvacuum_j1_leading_term(tau, z))
    assert relative_error_with_reference(j1, j1_scale * j1_series, j1) < 1e-4

    j2 = sl2_boundary_character(3, 2, tau, z, num_terms=500)
    j2_series = sl2_boundary_nonvacuum_j2_series_first_4_terms(tau, z)
    j2_scale = CC(j2) / CC(sl2_boundary_nonvacuum_j2_leading_term(tau, z))
    assert relative_error_with_reference(j2, j2_scale * j2_series, j2) < 3e-4

    sl3 = sl3_boundary_vacuum_character(tau, z1, z2, num_terms=500)
    assert relative_error(sl3, sl3_vacuum_series_first_2_terms(tau, z1, z2)) < 1e-4
