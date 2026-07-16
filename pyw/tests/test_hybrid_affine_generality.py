import pytest


@pytest.mark.sage
def test_b2_root_to_weight_uses_cartan_columns_and_round_trips():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra

    algebra = AffineLieAlgebra(["B", 2, 1])
    finite_root_lattice = algebra._finite_root_system.root_lattice()
    roots = tuple(finite_root_lattice.simple_roots().values()) + (
        finite_root_lattice.highest_root(),
    )

    for root in roots:
        weight = algebra.root_to_weight(root, finite=True)
        assert algebra.weight_to_root(weight, finite=True) == root
        assert weight.to_ambient() == root.to_ambient()


@pytest.mark.sage
def test_b2_fundamental_weight_form_is_symmetric_with_correct_root_lengths():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra

    algebra = AffineLieAlgebra(["B", 2, 1])
    fundamental_weights = algebra._finite_root_system.weight_space().fundamental_weights()
    simple_roots = algebra._finite_root_system.root_lattice().simple_roots()

    for first in fundamental_weights.values():
        for second in fundamental_weights.values():
            assert algebra.scalar_product(first, second) == algebra.scalar_product(second, first)
    assert algebra.scalar_product(simple_roots[1], simple_roots[1]) == 2
    assert algebra.scalar_product(simple_roots[2], simple_roots[2]) == 1


@pytest.mark.sage
@pytest.mark.parametrize(
    ("cartan_type", "expected_depth_one"),
    [(["A", 1, 1], 3), (["A", 2, 1], 7), (["B", 2, 1], 9)],
)
def test_vacuum_depth_one_is_finite_root_system_plus_zero(
    cartan_type,
    expected_depth_one,
):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import (
        build_bounded_weight_tree,
        classically_dominant_affine_weight,
    )

    algebra = AffineLieAlgebra(cartan_type)
    highest_weight = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    tree = build_bounded_weight_tree(highest_weight, order=1)
    weights = tree.weights_at_depth(1)

    assert len(weights) == expected_depth_one
    assert all(classically_dominant_affine_weight(weight) in weights for weight in weights)


@pytest.mark.sage
def test_weight_tree_rejects_nonintegral_finite_highest_weight():
    from sage.all import QQ

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.hybrid_affine_character import build_bounded_weight_tree

    algebra = AffineLieAlgebra(["A", 2, 1])
    highest_weight = QQ(1) / 2 * AffineWeight.affine_fundamental_weight(algebra, 1)

    with pytest.raises(NotImplementedError, match="dominant integral"):
        build_bounded_weight_tree(highest_weight, order=1)
