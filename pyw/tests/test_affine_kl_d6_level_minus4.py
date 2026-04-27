import pytest


@pytest.mark.sage
def test_d6_minus_4w0_task_prd_fixture_is_stable():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["D", 6, 1])
    fw = ala.fundamental_weights()
    lam = -4 * fw[0]
    kl_char = KazhdanLusztigCharacter(ala)

    Lambda_hat, w_to_Lambda = kl_char._find_dominant_Lambda(lam)

    assert Lambda_hat.to_legacy_expression() == "-Lambda[2] - Lambda[4] + 3*delta"
    assert tuple(int(i) for i in w_to_Lambda.reduced_word()) == (3, 1, 2, 0)
    assert ala.dim == 66


@pytest.mark.sage
def test_d6_minus_4w0_translation_enumerator_matches_bnb_for_task_order_zero():
    from sage.all import QQ

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["D", 6, 1])
    fw = ala.fundamental_weights()
    lam = -4 * fw[0]
    kl_char = KazhdanLusztigCharacter(ala)
    rho = ala.affine_rho()

    Lambda_hat, _ = kl_char._find_dominant_Lambda(lam)
    dominant = Lambda_hat + rho
    translation_order = 0 + max(0, int(QQ(dominant.grade) - QQ(lam.grade)))
    max_neg_shift = QQ(0) + QQ(dominant.grade)

    grid = kl_char._translations_by_n_shift(
        dominant,
        order=translation_order,
        max_neg_shift=max_neg_shift,
    )
    bnb = kl_char._translations_by_n_shift_bnb(
        dominant,
        order=translation_order,
        max_neg_shift=max_neg_shift,
    )

    def coeff_rows(translations):
        rows = []
        for t in translations:
            coeffs = t.translation_vector.monomial_coefficients()
            rows.append(
                tuple(
                    int(coeffs.get(i, 0))
                    for i in sorted(t.translation_vector.parent().index_set())
                )
            )
        return sorted(rows)

    assert len(grid) == 12
    assert coeff_rows(grid) == coeff_rows(bnb)
