import pytest


@pytest.mark.sage
def test_extended_affine_weyl_left_action_composition():
    """(g1*g2).action(x) == g1.action(g2.action(x))"""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import ExtendedAffineWeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])
    W = ExtendedAffineWeylGroup(ala)

    s0 = W.simple_reflection(0)
    s1 = W.simple_reflection(1)
    x = ala.affine_fundamental_weights()[1]

    assert (s0 * s1).action(x) == s0.action(s1.action(x))


@pytest.mark.sage
def test_extended_affine_weyl_s0_decomposition_matches_action():
    """s0 acts like s_theta * t_{-theta^vee} on AffineWeight."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import ExtendedAffineWeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])
    W = ExtendedAffineWeylGroup(ala)

    s0 = W.simple_reflection(0)
    x = ala.affine_fundamental_weights()[0]

    # manual: apply translation then finite reflection
    t = W.translation(-W.theta_coroot())
    manual = W._s_theta_weight.action(t.action(x).finite_part)
    expected = x.__class__(x.algebra, manual, x.level, t.action(x).grade)

    assert s0.action(x) == expected


@pytest.mark.sage
def test_translation_grade_uses_weight_coroot_pairing():
    """Translation should shift grade by -( <λ,β> + k|β|^2/2 )."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import ExtendedAffineWeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])
    W = ExtendedAffineWeylGroup(ala)

    x = ala.affine_fundamental_weights()[1]
    beta = -W._finite_coroot_space.simple_roots()[2]
    t = W.translation(beta)

    y = t.action(x)

    # expected pairing and norm are computed in ambient space
    pairing = x.finite_part.to_ambient().inner_product(beta.to_ambient())
    norm_sq = beta.norm_squared()
    expected_grade = x.grade - pairing - x.level * norm_sq / 2

    assert y.grade == expected_grade


@pytest.mark.sage
def test_affine_lie_algebra_affine_weyl_group_glue():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import ExtendedAffineWeylGroup, FiniteWeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])

    W_hat = ala.affine_weyl_group()
    assert isinstance(W_hat, ExtendedAffineWeylGroup)
    # cached
    assert ala.affine_weyl_group() is W_hat

    W = ala.finite_weyl_group()
    assert isinstance(W, FiniteWeylGroup)
    assert ala.finite_weyl_group() is W
