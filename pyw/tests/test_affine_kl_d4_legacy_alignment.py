import pytest


@pytest.mark.sage
def test_d4_minus_2w0_context_matches_legacy_getlambda_and_bruhat_flow():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["D", 4, 1])
    fw = ala.fundamental_weights()
    lam = -2 * fw[0]
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.build_context(lam, order=0)

    assert context.Lambda_hat.to_legacy_expression() == "-Lambda[2] + delta"
    assert tuple(int(i) for i in context.w_to_Lambda.reduced_word_list()) == (0,)
    assert tuple(int(i) for i in context.w_to_lambda.reduced_word_list()) == (0,)

    weyl_to_be_summed = context.weyl_to_be_summed()

    assert weyl_to_be_summed
    assert all(hasattr(w, "parent") for w in weyl_to_be_summed)


@pytest.mark.sage
def test_d4_minus_2w0_weyl_to_be_summed_avoids_affine_word_key_lookup():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["D", 4, 1])
    fw = ala.fundamental_weights()
    lam = -2 * fw[0]
    kl_char = KazhdanLusztigCharacter(ala)

    context = kl_char.build_context(lam, order=0)
    weyl_to_be_summed = context.weyl_to_be_summed()

    assert weyl_to_be_summed
    assert all(hasattr(w, "reduced_word") for w in weyl_to_be_summed)


@pytest.mark.sage
def test_d6_minus_4w0_fast_dominant_lambda_matches_reflection_chain_example():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["D", 6, 1])
    fw = ala.fundamental_weights()
    lam = -4 * fw[0]
    kl_char = KazhdanLusztigCharacter(ala)
    rho = ala.affine_rho()

    Lambda_hat, w_to_Lambda = kl_char._find_dominant_Lambda(lam)

    assert Lambda_hat.to_legacy_expression() == "-Lambda[2] - Lambda[4] + 3*delta"
    assert Lambda_hat.is_dominant() is False
    assert (Lambda_hat + rho).is_dominant()
    assert tuple(int(i) for i in w_to_Lambda.reduced_word_list()) == (3, 1, 2, 0)
