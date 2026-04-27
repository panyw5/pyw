import pytest


@pytest.mark.sage
def test_integrable_module_character_accepts_affine_weight_input():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)

    rep = IntegrableModuleCharacter(lam)

    assert rep.highest_weight() == lam.to_sagemath(extended=False)


@pytest.mark.sage
def test_affine_weight_to_sagemath_supports_affine_action_domain():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0) - AffineWeight.delta(ala)

    ordinary = lam.to_sagemath(extended=False)
    domain_weight = lam.to_sagemath()

    assert tuple(ordinary.to_vector()) == (1, 0, 0)
    assert tuple(domain_weight.to_vector()) == (1, 0, 0, -1)


@pytest.mark.sage
def test_integrable_module_character_strings_and_maximal_weights():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    rep = IntegrableModuleCharacter(lam)

    dmax = rep.dominant_maximal_weights()
    strings = rep.strings(3)

    assert dmax
    assert all(weight in strings for weight in dmax)
    assert all(len(values) >= 1 for values in strings.values())


@pytest.mark.sage
def test_integrable_module_character_depth_must_be_positive():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    rep = IntegrableModuleCharacter(lam)

    with pytest.raises(ValueError, match="positive integer"):
        rep.strings(0)


@pytest.mark.sage
def test_integrable_module_character_returns_truncated_q_series():
    from sage.all import var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    rep = IntegrableModuleCharacter(lam)

    q = var("q")
    character = rep.character(1)

    assert character != 0
    assert character.coefficient(q, 0) != 0
    assert character.coefficient(q, 2) == 0


@pytest.mark.sage
def test_integrable_module_character_preserves_a2_q_grading_regression():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    rep = IntegrableModuleCharacter(ala.fundamental_weights()[0])

    q = var("q")
    z1 = var("z1")
    z2 = var("z2")

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

    character = rep.character(3)

    assert character == expected
    assert character.coefficient(q, 0) == 1
    assert character.coefficient(q, 1) != 0
    assert character.coefficient(q, 2) != 0
    assert character.coefficient(q, 3) != 0


@pytest.mark.sage
def test_integrable_module_character_preserves_a2_strings_prefix_to_depth_six():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    rep = IntegrableModuleCharacter(ala.fundamental_weights()[0])

    strings = rep.strings(6)

    assert strings == {rep.highest_weight(): [1, 2, 5, 10, 20, 36]}


@pytest.mark.sage
def test_integrable_module_character_preserves_a2_q_grading_to_order_five_regression():
    from sage.all import SR, var

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    rep = IntegrableModuleCharacter(ala.fundamental_weights()[0])

    q = var("q")
    z1 = var("z1")
    z2 = var("z2")

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
    expected += (
        z1**2 * z2**2
        + 2 * z1**3
        + 2 * z2**3
        + z1**4 / z2**2
        + 10 * z1 * z2
        + z2**4 / z1**2
        + 10 * z1**2 / z2
        + 10 * z2**2 / z1
        + 2 * z1**3 / z2**3
        + 2 * z2**3 / z1**3
        + 10 * z1 / z2**2
        + 10 * z2 / z1**2
        + z1**2 / z2**4
        + 10 / (z1 * z2)
        + z2**2 / z1**4
        + 2 / z1**3
        + 2 / z2**3
        + 1 / (z1**2 * z2**2)
        + 20
    ) * q**4
    expected += (
        2 * z1**2 * z2**2
        + 5 * z1**3
        + 5 * z2**3
        + 2 * z1**4 / z2**2
        + 20 * z1 * z2
        + 2 * z2**4 / z1**2
        + 20 * z1**2 / z2
        + 20 * z2**2 / z1
        + 5 * z1**3 / z2**3
        + 5 * z2**3 / z1**3
        + 20 * z1 / z2**2
        + 20 * z2 / z1**2
        + 2 * z1**2 / z2**4
        + 20 / (z1 * z2)
        + 2 * z2**2 / z1**4
        + 5 / z1**3
        + 5 / z2**3
        + 2 / (z1**2 * z2**2)
        + 36
    ) * q**5

    character = rep.character(5)

    assert character == expected
    assert character.coefficient(q, 5) == expected.coefficient(q, 5)


@pytest.mark.sage
def test_integrable_module_character_accepts_explicit_translations():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    rep = IntegrableModuleCharacter(lam)

    explicit = [0]
    character = rep.character(0, manual_translations=explicit)

    assert character != 0


@pytest.mark.sage
def test_integrable_module_character_auto_path_does_not_instantiate_kl_character(monkeypatch):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter, KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    rep = IntegrableModuleCharacter(lam)

    def _forbidden_init(self, algebra):
        raise AssertionError("IntegrableModuleCharacter should not instantiate KazhdanLusztigCharacter")

    monkeypatch.setattr(KazhdanLusztigCharacter, "__init__", _forbidden_init)

    character = rep.character(0)

    assert character != 0
