import pytest


def _d4_vacuum():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight

    algebra = AffineLieAlgebra(["D", 4, 1])
    highest_weight = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    return algebra, highest_weight


@pytest.mark.sage
def test_bounded_d4_vacuum_weight_tree_through_depth_two():
    from pyw.core.hybrid_affine_character import build_bounded_weight_tree

    _, highest_weight = _d4_vacuum()
    tree = build_bounded_weight_tree(highest_weight, order=2)

    assert [len(tree.weights_at_depth(depth)) for depth in range(3)] == [1, 25, 169]
    assert all(weight.grade == -depth for depth in range(3) for weight in tree.weights_at_depth(depth))


@pytest.mark.sage
def test_bounded_weight_tree_accepts_general_finite_dominant_highest_weight():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import build_bounded_weight_tree

    algebra = AffineLieAlgebra(["A", 1, 1])
    nonvacuum = AffineWeight.affine_fundamental_weight(algebra, 1)
    shifted = nonvacuum + 3 * algebra.affine_delta()

    tree = build_bounded_weight_tree(shifted, order=1)

    assert len(tree.weights_at_depth(0)) == 2
    assert all(weight.grade == 3 for weight in tree.weights_at_depth(0))
    assert tree.weights_at_depth(1)
    assert all(weight.grade == 2 for weight in tree.weights_at_depth(1))


@pytest.mark.sage
def test_freudenthal_denominator_detects_known_d4_targets():
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import freudenthal_denominator

    algebra, highest_weight = _d4_vacuum()
    finite_weights = algebra._finite_root_system.weight_space().fundamental_weights()
    nondegenerate = AffineWeight(algebra, finite_weights[2], -2, -1)
    degenerate = AffineWeight(algebra, 2 * finite_weights[1], -2, -2)

    assert freudenthal_denominator(highest_weight, highest_weight) == 0
    assert freudenthal_denominator(highest_weight, nondegenerate) != 0
    assert freudenthal_denominator(highest_weight, degenerate) == 0


@pytest.mark.sage
def test_freudenthal_routes_degenerate_targets_and_matches_low_d4_multiplicities():
    from pyw.core.affine_weight import affine_weight_key
    from pyw.core.hybrid_affine_character import (
        build_bounded_weight_tree,
        compute_bounded_freudenthal_multiplicities,
        freudenthal_denominator,
    )

    _, highest_weight = _d4_vacuum()
    tree = build_bounded_weight_tree(highest_weight, order=2)
    degenerate_targets = []

    def degenerate_multiplicity(weight):
        degenerate_targets.append(weight)
        return 2

    multiplicities = compute_bounded_freudenthal_multiplicities(
        tree,
        degenerate_multiplicity=degenerate_multiplicity,
    )
    by_finite_labels = {
        (
            depth,
            tuple(weight.dynkin_labels()[index] for index in range(1, 5)),
        ): multiplicities[affine_weight_key(weight)]
        for depth in range(3)
        for weight in tree.weights_at_depth(depth)
        if affine_weight_key(weight) in multiplicities
    }

    assert len(degenerate_targets) == 3
    assert all(
        freudenthal_denominator(highest_weight, weight) == 0
        for weight in degenerate_targets
    )
    assert by_finite_labels[(1, (0, 1, 0, 0))] == 1
    assert by_finite_labels[(2, (0, 2, 0, 0))] == 1
    assert by_finite_labels[(2, (0, 1, 0, 0))] == 6
    assert by_finite_labels[(2, (0, 0, 0, 0))] == 17


@pytest.mark.sage
def test_freudenthal_rejects_nonintegral_degenerate_result():
    from sage.all import QQ

    from pyw.core.hybrid_affine_character import (
        build_bounded_weight_tree,
        compute_bounded_freudenthal_multiplicities,
    )

    _, highest_weight = _d4_vacuum()
    tree = build_bounded_weight_tree(highest_weight, order=2)

    with pytest.raises(ArithmeticError, match="nonnegative integer"):
        compute_bounded_freudenthal_multiplicities(
            tree,
            degenerate_multiplicity=lambda _weight: QQ(1) / 2,
        )


@pytest.mark.sage
@pytest.mark.slow
def test_d4_order_four_freudenthal_matches_character_dimensions():
    from collections import Counter

    from pyw.core.affine_weight import affine_weight_key
    from pyw.core.hybrid_affine_character import (
        build_bounded_weight_tree,
        classically_dominant_affine_weight,
        compute_bounded_freudenthal_multiplicities,
    )

    _, highest_weight = _d4_vacuum()
    tree = build_bounded_weight_tree(highest_weight, order=4)
    expected_degenerate_values = {
        (2, (2, 0, 0, 0)): 2,
        (2, (0, 0, 2, 0)): 2,
        (2, (0, 0, 0, 2)): 2,
        (3, (1, 0, 1, 1)): 8,
        (4, (2, 1, 0, 0)): 16,
        (4, (0, 1, 2, 0)): 16,
        (4, (0, 1, 0, 2)): 16,
    }
    degenerate_calls = []

    def degenerate_multiplicity(weight):
        key = (
            -int(weight.grade),
            tuple(int(weight.dynkin_labels()[index]) for index in range(1, 5)),
        )
        degenerate_calls.append(key)
        return expected_degenerate_values[key]

    multiplicities = compute_bounded_freudenthal_multiplicities(
        tree,
        degenerate_multiplicity=degenerate_multiplicity,
    )
    dimensions = [
        sum(
            multiplicities[
                affine_weight_key(classically_dominant_affine_weight(weight))
            ]
            for weight in tree.weights_at_depth(depth)
        )
        for depth in range(5)
    ]
    dominant_counts = Counter(
        -int(tree.weight_by_key[key].grade) for key in multiplicities
    )

    assert [len(tree.weights_at_depth(depth)) for depth in range(5)] == [
        1,
        25,
        169,
        625,
        1681,
    ]
    assert dominant_counts == {0: 1, 1: 2, 2: 7, 3: 15, 4: 30}
    assert set(degenerate_calls) == set(expected_degenerate_values)
    assert len(degenerate_calls) == 7
    assert dimensions == [1, 28, 329, 2632, 16380]
