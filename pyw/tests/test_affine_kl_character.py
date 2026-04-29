import pytest


@pytest.mark.sage
def test_kl_character_builds_context():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.prepare_data(lam, order=1)

    assert context.Lambda_hat.is_dominant()
    assert context.quotient_representatives


@pytest.mark.sage
def test_kl_character_numerator_terms_are_available():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    terms = kl_char.numerator_terms(lam, order=1)

    assert terms
    assert all(term.coefficient != 0 for term in terms)


@pytest.mark.sage
def test_kl_character_numerator_terms_accept_explicit_translations():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)
    beta = ala.affine_weyl_group()._finite_coroot_space.simple_roots()[1]

    terms = kl_char.numerator_terms(lam, order=1, manual_translations=[0, beta])

    assert terms
    assert all(term.coefficient != 0 for term in terms)


@pytest.mark.sage
def test_kl_character_weyl_to_be_summed_matches_direct_bruhat_filter():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.bruhat import BruhatOrder
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.prepare_data(lam, order=1)
    bruhat = BruhatOrder(context.W_affine_as_words[0].parent())
    direct = [w for w in context.quotient_representatives if bruhat.le(context.w_to_lambda, w)]
    direct_words = sorted(tuple(int(i) for i in w.reduced_word()) for w in direct)

    optimized = context.weyl_to_be_summed()
    optimized_words = sorted(tuple(int(i) for i in w.reduced_word()) for w in optimized)

    assert optimized_words == direct_words


@pytest.mark.sage
def test_character_numerator_legacy_matches_numerator_terms_a21():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    terms = kl_char.numerator_terms(lam, order=1)
    legacy = kl_char.character_numerator_legacy(lam, order=1)

    terms_map = {term.weight: term.coefficient for term in terms}
    legacy_map = {}
    for entry in legacy:
        for weight, coeff in entry.items():
            legacy_map[weight] = coeff

    assert len(terms_map) == len(legacy_map), (
        f"term count mismatch: numerator_terms={len(terms_map)}, "
        f"character_numerator_legacy={len(legacy_map)}"
    )

    for weight, coeff in terms_map.items():
        assert weight in legacy_map, f"weight {weight} missing from character_numerator_legacy"
        assert legacy_map[weight] == coeff, (
            f"coefficient mismatch for weight {weight}: "
            f"numerator_terms={coeff}, character_numerator_legacy={legacy_map[weight]}"
        )
