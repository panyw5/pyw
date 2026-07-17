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
    int_char = IntegrableModuleCharacter(ala)

    q = var("q")
    rename = {var("z1"): var("b1"), var("z2"): var("b2")}

    for order in range(4):
        kl_result = kl_char.character(lam, order=order)
        int_result = int_char.character(lam, order).subs(rename)

        for grade in range(order + 1):
            kl_coeff = kl_result.coefficient(q, grade)
            int_coeff = int_result.coefficient(q, grade)
            assert (kl_coeff - int_coeff).expand() == 0


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
    b1 = var("b1")
    b2 = var("b2")

    result = kl_char.character(lam, order=3)

    expected = SR(1)
    expected += (
        b1 * b2 + b1**2 / b2 + b2**2 / b1 + b1 / b2**2 + b2 / b1**2 + 1 / (b1 * b2) + 2
    ) * q
    expected += (
        2 * b1 * b2
        + 2 * b1**2 / b2
        + 2 * b2**2 / b1
        + 2 * b1 / b2**2
        + 2 * b2 / b1**2
        + 2 / (b1 * b2)
        + 5
    ) * q**2
    expected += (
        b1**3
        + b2**3
        + 5 * b1 * b2
        + 5 * b1**2 / b2
        + 5 * b2**2 / b1
        + b1**3 / b2**3
        + b2**3 / b1**3
        + 5 * b1 / b2**2
        + 5 * b2 / b1**2
        + 5 / (b1 * b2)
        + 1 / b1**3
        + 1 / b2**3
        + 10
    ) * q**3

    assert (result - expected).expand() == 0


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
def test_kl_character_accepts_show_progress_and_debug():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    kl_char = KazhdanLusztigCharacter(ala)

    result_default = kl_char.character(lam, order=1)
    result_progress = kl_char.character(lam, order=1, show_progress=True)
    result_debug = kl_char.character(lam, order=1, debug=True)
    result_both = kl_char.character(lam, order=1, show_progress=True, debug=True)

    assert result_default == result_progress
    assert result_default == result_debug
    assert result_default == result_both
