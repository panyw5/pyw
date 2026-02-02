import pytest


# Test multiple Cartan types to ensure equations hold universally
AFFINE_CARTAN_TYPES = [
    ["A", 2, 1],  # A2^(1) - simply laced
    ["A", 3, 1],  # A3^(1) - simply laced, higher rank
    ["B", 2, 1],  # B2^(1) - non-simply laced
    ["C", 2, 1],  # C2^(1) - non-simply laced
    ["D", 4, 1],  # D4^(1) - simply laced with triality
    ["G", 2, 1],  # G2^(1) - exceptional, non-simply laced
    ["F", 4, 1],  # F4^(1) - exceptional
    ["E", 6, 1],  # E6^(1) - exceptional, simply laced
]


@pytest.mark.sage
@pytest.mark.parametrize("cartan_type", AFFINE_CARTAN_TYPES)
def test_extended_affine_weyl_left_action_composition(cartan_type):
    """(g1*g2).action(x) == g1.action(g2.action(x)) for all affine types."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import AffineWeylGroupSemidirect

    ala = AffineLieAlgebra(cartan_type)
    W = AffineWeylGroupSemidirect(ala)

    # Test with s0 and s1
    s0 = W.simple_reflection(0)
    s1 = W.simple_reflection(1)
    x = ala.affine_fundamental_weights()[1]

    assert (s0 * s1).action(x) == s0.action(s1.action(x))

    # Test with other combinations if rank >= 2
    index_set = list(ala.root_system().index_set())
    if len([i for i in index_set if i != 0]) >= 2:
        finite_indices = [i for i in index_set if i != 0]
        i1, i2 = finite_indices[0], finite_indices[1]
        si1 = W.simple_reflection(i1)
        si2 = W.simple_reflection(i2)

        assert (si1 * si2).action(x) == si1.action(si2.action(x))


@pytest.mark.sage
@pytest.mark.parametrize("cartan_type", AFFINE_CARTAN_TYPES)
def test_extended_affine_weyl_s0_decomposition_matches_action(cartan_type):
    """s0 = s_theta * t_{-theta^vee} for all affine types."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import AffineWeylGroupSemidirect

    ala = AffineLieAlgebra(cartan_type)
    W = AffineWeylGroupSemidirect(ala)

    s0 = W.simple_reflection(0)
    x = ala.affine_fundamental_weights()[0]

    # manual: apply translation then finite reflection
    t = W.translation(-W.theta_coroot())
    manual = W._s_theta_weight.action(t.action(x).finite_part)
    expected = x.__class__(x.algebra, manual, x.level, t.action(x).grade)

    assert s0.action(x) == expected


@pytest.mark.sage
@pytest.mark.parametrize("cartan_type", AFFINE_CARTAN_TYPES)
def test_translation_grade_uses_weight_coroot_pairing(cartan_type):
    """Translation shifts grade by -(<λ,β> + k|β|^2/2) for all affine types."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import AffineWeylGroupSemidirect

    ala = AffineLieAlgebra(cartan_type)
    W = AffineWeylGroupSemidirect(ala)

    # Test with first non-zero fundamental weight
    index_set = list(ala.root_system().index_set())
    finite_indices = [i for i in index_set if i != 0]
    test_index = finite_indices[0] if finite_indices else 1

    x = ala.affine_fundamental_weights()[test_index]

    # Test with a simple coroot
    simple_coroots = W._finite_coroot_space.simple_roots()
    beta = simple_coroots[test_index]
    t = W.translation(beta)

    y = t.action(x)

    # expected pairing and norm are computed in ambient space
    pairing = x.finite_part.to_ambient().inner_product(beta.to_ambient())
    norm_sq = beta.norm_squared()
    expected_grade = x.grade - pairing - x.level * norm_sq / 2

    assert y.grade == expected_grade

    # Test with negative coroot as well
    beta_neg = -beta
    t_neg = W.translation(beta_neg)
    y_neg = t_neg.action(x)

    pairing_neg = x.finite_part.to_ambient().inner_product(beta_neg.to_ambient())
    norm_sq_neg = beta_neg.norm_squared()
    expected_grade_neg = x.grade - pairing_neg - x.level * norm_sq_neg / 2

    assert y_neg.grade == expected_grade_neg


@pytest.mark.sage
def test_affine_lie_algebra_affine_weyl_group_glue():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.weyl_group import AffineWeylGroupSemidirect, FiniteWeylGroup

    ala = AffineLieAlgebra(["A", 2, 1])

    W_hat = ala.affine_weyl_group()
    assert isinstance(W_hat, AffineWeylGroupSemidirect)
    # cached
    assert ala.affine_weyl_group() is W_hat

    W = ala.finite_weyl_group()
    assert isinstance(W, FiniteWeylGroup)
    assert ala.finite_weyl_group() is W
