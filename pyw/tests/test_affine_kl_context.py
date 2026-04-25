import pytest
from sage.all import QQ


def _coeff_tuple(ala, beta):
    idxs = [int(i) for i in ala.affine_weyl_group()._finite_coroot_space.index_set()]
    coeffs = beta.monomial_coefficients()
    return tuple(int(coeffs.get(i, 0)) for i in idxs)


def _neg_shift_formula(weight, beta):
    lam_amb = weight.finite_part.to_ambient()
    beta_amb = beta.to_ambient()
    pairing = QQ(lam_amb.inner_product(beta_amb))
    beta_norm_sq = QQ(beta_amb.inner_product(beta_amb))
    return pairing + QQ(weight.level) * beta_norm_sq / QQ(2)


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

    context = kl_char.prepare_data(lambda_hat, order=1)

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

    context = kl_char.prepare_data(lambda_hat, order=1)

    assert any(tuple(int(i) for i in w.reduced_word()) == () for w in context.stabilizer_candidates)


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

    context = kl_char.prepare_data(lambda_hat, order=1)

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

    context = kl_char.prepare_data(lambda_hat, order=2)

    assert context.translation_bounds is None
    assert len(context.translations) == 3
    assert context.translations[0].translation_vector == ala.affine_weyl_group()._zero_beta
    assert any(t.translation_vector == -beta for t in context.translations)
    assert any(
        t.translation_vector == -(beta + ala.affine_weyl_group()._finite_coroot_space.simple_roots()[2])
        for t in context.translations
    )


@pytest.mark.sage
def test_translations_by_n_shift_returns_translation_operators():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    translations = kl_char._translations_by_n_shift(ala.affine_rho(), order=1)

    assert translations
    assert all(hasattr(t, "translation_vector") for t in translations)
    assert all(t.finite_part.reduced_word() == [] for t in translations)


@pytest.mark.sage
def test_translations_by_n_shift_matches_direct_inequality_filter_in_bounded_box():
    from itertools import product

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    weight = ala.affine_rho()
    order = 30
    method = kl_char._translations_by_n_shift(weight, order=order)

    semidirect = ala.affine_weyl_group()
    idxs = [int(i) for i in semidirect._finite_coroot_space.index_set()]
    basis = semidirect._finite_coroot_space.simple_roots()
    box = range(-8, 9)

    oracle = set()
    for coeffs in product(box, repeat=len(idxs)):
        beta = semidirect._zero_beta
        for i, c in zip(idxs, coeffs):
            if c:
                beta += int(c) * basis[i]
        neg_shift = _neg_shift_formula(weight, beta)
        if QQ(0) <= neg_shift <= QQ(order) + QQ(weight.grade):
            oracle.add(tuple(int(c) for c in coeffs))

    method_coeffs = {_coeff_tuple(ala, t.translation_vector) for t in method}
    assert oracle.issubset(method_coeffs)


@pytest.mark.sage
def test_translations_by_n_shift_finds_coefficients_beyond_legacy_box():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    weight = ala.affine_rho()
    translations = kl_char._translations_by_n_shift(weight, order=200)

    coeffs = [_coeff_tuple(ala, t.translation_vector) for t in translations]
    assert any(max(abs(c) for c in row) > 5 for row in coeffs)


@pytest.mark.sage
def test_translations_by_n_shift_bnb_matches_reference_enumerator():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    weight = ala.affine_rho()

    ref = kl_char._translations_by_n_shift(weight, order=40)
    alt = kl_char._translations_by_n_shift_bnb(weight, order=40)

    ref_coeffs = {_coeff_tuple(ala, t.translation_vector) for t in ref}
    alt_coeffs = {_coeff_tuple(ala, t.translation_vector) for t in alt}
    assert alt_coeffs == ref_coeffs


@pytest.mark.sage
def test_translations_by_n_shift_bnb_reports_pruning_stats():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    weight = ala.affine_rho()

    result = kl_char._translations_by_n_shift_bnb(weight, order=80, return_stats=True)
    stats = result["stats"]

    assert stats["box_points"] > 0
    assert stats["accepted"] > 0
    assert stats["visited_leaves"] <= stats["box_points"]
    assert stats["pruned_branches"] > 0


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

    context = kl_char.prepare_data(
        lambda_hat,
        order=1,
        translations=[W.translation(beta), beta, 0, beta],
    )

    assert context.translation_bounds is None
    assert [t.translation_vector for t in context.translations] == [W._zero_beta, beta]
