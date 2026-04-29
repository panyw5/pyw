import pytest


@pytest.mark.sage
def test_character_numerator_legacy_a21():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    legacy = kl_char.character_numerator_legacy(lam, order=1)

    assert legacy
    assert all(len(entry) == 1 for entry in legacy)
    assert all(list(entry.values())[0] != 0 for entry in legacy)


@pytest.mark.sage
def test_character_numerator_legacy_accepts_explicit_translations():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)
    beta = ala.affine_weyl_group()._finite_coroot_space.simple_roots()[1]

    legacy = kl_char.character_numerator_legacy(lam, order=1, translations=[0, beta])

    assert legacy
    assert all(list(entry.values())[0] != 0 for entry in legacy)


@pytest.mark.sage
def test_character_a21():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    result = kl_char.character(lam, order=1)

    assert result != 0
