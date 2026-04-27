import pytest


@pytest.mark.sage
def test_integrable_module_character_accepts_affine_weight_input():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)

    rep = IntegrableModuleCharacter(lam)

    assert rep.highest_weight() == lam.to_sagemath()


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
def test_integrable_module_character_accepts_explicit_translations():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import IntegrableModuleCharacter

    ala = AffineLieAlgebra(["A", 1, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 0)
    rep = IntegrableModuleCharacter(lam)

    explicit = [0]
    character = rep.character(0, translations=explicit)

    assert character != 0
