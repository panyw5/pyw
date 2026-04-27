import pytest


@pytest.mark.sage
def test_kl_character_returns_symbolic_expression():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    kl_char = KazhdanLusztigCharacter(ala)

    result = kl_char.character(lam, order=1)

    assert result in SR
    q = var("q")
    assert result.coefficient(q, 0) == 1


@pytest.mark.sage
def test_kl_character_matches_integrable_module_character():
    from sage.all import var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter, KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)

    kl_char = KazhdanLusztigCharacter(ala)
    int_char = IntegrableModuleCharacter(lam)

    q = var("q")

    for order in range(4):
        kl_result = kl_char.character(lam, order=order)
        int_result = int_char.character(order)

        for grade in range(order + 1):
            kl_coeff = kl_result.coefficient(q, grade)
            int_coeff = int_result.coefficient(q, grade)
            assert kl_coeff == int_coeff


@pytest.mark.sage
def test_kl_character_preserves_a2_q_grading():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    kl_char = KazhdanLusztigCharacter(ala)

    q = var("q")
    z1 = var("z1")
    z2 = var("z2")

    result = kl_char.character(lam, order=3)

    expected = SR(1)
    expected += (z1 * z2 + z1**2 / z2 + z2**2 / z1 + z1 / z2**2 + z2 / z1**2 + 1 / (z1 * z2) + 2) * q
    expected += (2 * z1 * z2 + 2 * z1**2 / z2 + 2 * z2**2 / z1 + 2 * z1 / z2**2 + 2 * z2 / z1**2 + 2 / (z1 * z2) + 5) * q**2
    expected += (
        z1**3
        + z2**3
        + 5 * z1 * z2
        + 5 * z1**2 / z2
        + 5 * z2**2 / z1
        + z1**3 / z2**3
        + z2**3 / z1**3
        + 5 * z1 / z2**2
        + 5 * z2 / z1**2
        + 5 / (z1 * z2)
        + 1 / z1**3
        + 1 / z2**3
        + 10
    ) * q**3

    assert result == expected


@pytest.mark.sage
def test_kl_character_accepts_manual_translations():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    kl_char = KazhdanLusztigCharacter(ala)

    result = kl_char.character(lam, order=0, manual_translations=[0])

    assert result != 0


@pytest.mark.sage
def test_kl_formal_character_differs_from_character():
    from sage.all import var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    kl_char = KazhdanLusztigCharacter(ala)

    q = var("q")
    z1 = var("z1")
    z2 = var("z2")

    formal_result = kl_char.formal_character(lam, order=3)
    ws_result = kl_char.character(lam, order=3)

    for grade in range(4):
        formal_coeff = formal_result[grade]
        ws_coeff = ws_result.coefficient(q, grade)

        if grade > 0:
            ws_dim = ws_coeff.subs({z1: 1, z2: 1})
            assert formal_coeff != ws_dim
