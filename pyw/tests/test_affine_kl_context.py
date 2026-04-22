import pytest


@pytest.mark.sage
def test_affine_kl_context_finds_dominant_conjugate_identity_case():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials
    from sage.all import WeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    _ = KazhdanLusztigPolynomials(WeylGroup(["A", 2, 1]))
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.build_context(lambda_hat, order=1)

    assert context.Lambda_hat.is_dominant()
    assert context.quotient_weight(context.w_to_Lambda) == context.Lambda_hat
    assert context.apply(context.w_to_lambda, context.Lambda_hat + context.rho_hat) - context.rho_hat == context.lambda_hat


@pytest.mark.sage
def test_affine_kl_context_stabilizer_contains_identity():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials
    from sage.all import WeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    _ = KazhdanLusztigPolynomials(WeylGroup(["A", 2, 1]))
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.build_context(lambda_hat, order=1)

    assert any(len(w.reduced_word()) == 0 for w in context.stabilizer_candidates)


@pytest.mark.sage
def test_affine_kl_context_quotient_representatives_are_unique_by_weight():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials
    from sage.all import WeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    _ = KazhdanLusztigPolynomials(WeylGroup(["A", 2, 1]))
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.build_context(lambda_hat, order=1)

    seen = set()
    for w in context.quotient_representatives:
        image = context.quotient_weight(w)
        key = tuple(sorted(image.dynkin_labels().items())) + ((-1, image.grade),)
        assert key not in seen
        seen.add(key)


@pytest.mark.sage
def test_affine_kl_context_defaults_to_exact_translation_selection():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)
    beta = ala.affine_weyl_group()._finite_coroot_space.simple_roots()[1]

    context = kl_char.build_context(lambda_hat, order=2)

    assert context.translation_bounds is None
    assert len(context.translations) == 3
    assert context.translations[0].translation_vector == ala.affine_weyl_group()._zero_beta
    assert any(t.translation_vector == -beta for t in context.translations)
    assert any(
        t.translation_vector == -(beta + ala.affine_weyl_group()._finite_coroot_space.simple_roots()[2])
        for t in context.translations
    )


@pytest.mark.sage
def test_exact_translations_by_n_shift_returns_translation_operators():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    translations = kl_char._exact_translations_by_n_shift(ala.affine_rho(), order=1)

    assert translations
    assert all(hasattr(t, "translation_vector") for t in translations)
    assert all(t.finite_part.reduced_word() == [] for t in translations)


@pytest.mark.sage
def test_affine_kl_context_translation_bounds_remain_available_as_override():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.build_context(
        lambda_hat,
        order=1,
        translation_bounds={1: (-1, 1), 2: (0, 0)},
    )

    W = ala.affine_weyl_group()
    assert context.translation_bounds == {1: (-1, 1), 2: (0, 0)}
    assert [t.translation_vector for t in context.translations] == [
        W._zero_beta,
        W._finite_coroot_space.simple_roots()[1] * -1,
        W._finite_coroot_space.simple_roots()[1],
    ]


@pytest.mark.sage
def test_affine_kl_context_rejects_conflicting_translation_inputs():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)
    beta = ala.affine_weyl_group()._finite_coroot_space.simple_roots()[1]

    with pytest.raises(ValueError, match="translation_bounds or translations"):
        kl_char.build_context(
            lambda_hat,
            order=1,
            translation_bounds={1: (0, 0)},
            translations=[beta],
        )


@pytest.mark.sage
def test_affine_kl_context_stores_normalized_translation_set():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)
    W = ala.affine_weyl_group()
    beta = W._finite_coroot_space.simple_roots()[1]

    context = kl_char.build_context(
        lambda_hat,
        order=1,
        translations=[W.translation(beta), beta, 0, beta],
    )

    assert context.translation_bounds is None
    assert [t.translation_vector for t in context.translations] == [W._zero_beta, beta]
