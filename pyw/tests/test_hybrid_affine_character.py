import pytest
from sage.all import QQ


@pytest.mark.sage
@pytest.mark.parametrize("cartan_type", [["A", 2, 1], ["D", 4, 1], ["B", 2, 1]])
def test_affine_simple_root_coordinates_round_trip(cartan_type):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import affine_simple_root_coordinates

    algebra = AffineLieAlgebra(cartan_type)
    simple_roots = algebra.affine_simple_roots()
    expected = tuple(QQ(index + 1) for index in range(algebra.rank + 1))
    difference = AffineWeight.zero(algebra)
    for index, coefficient in enumerate(expected):
        difference += coefficient * simple_roots[index]

    assert affine_simple_root_coordinates(algebra, difference) == expected


@pytest.mark.sage
def test_affine_positive_root_cone_boundaries():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import is_in_affine_positive_root_cone

    algebra = AffineLieAlgebra(["A", 2, 1])
    simple_roots = algebra.affine_simple_roots()

    assert is_in_affine_positive_root_cone(algebra, AffineWeight.zero(algebra))
    assert is_in_affine_positive_root_cone(algebra, simple_roots[0] + simple_roots[2])
    assert not is_in_affine_positive_root_cone(algebra, -simple_roots[1])
    assert not is_in_affine_positive_root_cone(algebra, QQ(1) / 2 * simple_roots[1])
    assert not is_in_affine_positive_root_cone(
        algebra, AffineWeight.affine_fundamental_weight(algebra, 0)
    )


@pytest.mark.sage
def test_affine_weight_key_is_exact_and_parent_independent():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight, affine_weight_key

    algebra = AffineLieAlgebra(["A", 2, 1])
    weight = AffineWeight.affine_fundamental_weight(algebra, 1)
    rebuilt = algebra.from_sagemath(weight.to_sagemath(extended=True))

    assert affine_weight_key(weight) == affine_weight_key(rebuilt)
    shifted = weight + QQ(1) / 10**20 * algebra.affine_delta()
    assert affine_weight_key(shifted) != affine_weight_key(weight)


@pytest.mark.sage
def test_non_simply_laced_affine_dynkin_labels_use_comarks():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight

    algebra = AffineLieAlgebra(["B", 2, 1])
    weight = AffineWeight.affine_fundamental_weight(algebra, 2)

    assert weight.dynkin_labels() == {0: 0, 1: 0, 2: 1}
    assert list(weight.to_sagemath(extended=False).to_vector()) == [0, 0, 1]


@pytest.mark.sage
def test_from_sagemath_reads_grade_only_from_extended_weights():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight

    algebra = AffineLieAlgebra(["A", 2, 1])
    weight = AffineWeight.affine_fundamental_weight(algebra, 2)

    ordinary = algebra.from_sagemath(weight.to_sagemath(extended=False))
    extended = algebra.from_sagemath(weight.to_sagemath(extended=True))

    assert ordinary.grade == 0
    assert extended.grade == weight.grade
    assert ordinary == weight
    assert extended == weight


@pytest.mark.sage
def test_bounded_kl_orbit_has_exact_bidirectional_indexes():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight, affine_weight_key
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)
    orbit = KazhdanLusztigCharacter(algebra).prepare_bounded_kl_orbit(lambda_hat, order=1)

    assert orbit.lambda_hat == lambda_hat
    assert isinstance(orbit.representatives, tuple)
    assert isinstance(orbit.stabilizer, tuple)
    for weyl_key, weight in orbit.weyl_to_weight.items():
        assert orbit.weight_to_weyl[affine_weight_key(weight)].reduced_word() == list(weyl_key)
        assert affine_weight_key(weight) in orbit.by_grade[weight.grade]

    with pytest.raises(TypeError):
        orbit.weyl_to_weight[()] = lambda_hat
    with pytest.raises(TypeError):
        orbit.weight_to_weyl[affine_weight_key(lambda_hat)] = orbit.representatives[0]


@pytest.mark.sage
def test_candidates_above_apply_affine_root_cone_filter():
    from types import MappingProxyType

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight, affine_weight_key
    from pyw.core.character import BoundedKLOrbit
    from pyw.core.weyl_group import element_word_list

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)
    simple_roots = algebra.affine_simple_roots()
    affine_weyl_group = algebra.affine_weyl_group_sage()
    representatives = (
        affine_weyl_group.one(),
        affine_weyl_group.from_reduced_word([0]),
        affine_weyl_group.from_reduced_word([1]),
    )
    orbit_weights = (
        lambda_hat,
        lambda_hat + simple_roots[0] + simple_roots[1],
        lambda_hat - simple_roots[0],
    )
    weyl_to_weight = {
        tuple(element_word_list(representative)): weight
        for representative, weight in zip(representatives, orbit_weights, strict=True)
    }
    orbit = BoundedKLOrbit(
        lambda_hat=lambda_hat,
        Lambda_hat=lambda_hat,
        w_to_lambda=affine_weyl_group.one(),
        stabilizer=(affine_weyl_group.one(),),
        representatives=representatives,
        weyl_to_weight=MappingProxyType(weyl_to_weight),
        weight_to_weyl=MappingProxyType(
            {
                affine_weight_key(weight): representative
                for representative, weight in zip(
                    representatives, orbit_weights, strict=True
                )
            }
        ),
        by_grade=MappingProxyType({}),
    )

    candidates = orbit.candidates_above(lambda_hat)

    assert [candidate.representative for candidate in candidates] == list(representatives[:2])
    assert [candidate.difference_coordinates for candidate in candidates] == [
        (0, 0, 0),
        (1, 1, 0),
    ]


@pytest.mark.sage
@pytest.mark.slow
def test_d4_order_four_candidate_support_benchmark():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["D", 4, 1])
    finite_fundamental_weights = algebra._finite_root_system.weight_space().fundamental_weights()
    lambda_hat = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    targets = (
        AffineWeight(algebra, 2 * finite_fundamental_weights[1], -2, -2),
        AffineWeight(algebra, 2 * finite_fundamental_weights[3], -2, -2),
        AffineWeight(algebra, 2 * finite_fundamental_weights[4], -2, -2),
        AffineWeight(
            algebra,
            finite_fundamental_weights[1]
            + finite_fundamental_weights[3]
            + finite_fundamental_weights[4],
            -2,
            -3,
        ),
        AffineWeight(
            algebra,
            2 * finite_fundamental_weights[1] + finite_fundamental_weights[2],
            -2,
            -4,
        ),
        AffineWeight(
            algebra,
            2 * finite_fundamental_weights[3] + finite_fundamental_weights[2],
            -2,
            -4,
        ),
        AffineWeight(
            algebra,
            2 * finite_fundamental_weights[4] + finite_fundamental_weights[2],
            -2,
            -4,
        ),
    )

    orbit = KazhdanLusztigCharacter(algebra).prepare_bounded_kl_orbit(
        lambda_hat,
        order=4,
    )
    supports = tuple(orbit.candidates_above(target) for target in targets)

    assert len(orbit.representatives) == 2112
    assert [len(support) for support in supports] == [12, 12, 12, 28, 58, 58, 58]
    assert len({candidate.weyl_key for support in supports for candidate in support}) == 81
