import pytest


@pytest.mark.sage
def test_bounded_affine_positive_roots_include_real_and_imaginary_roots():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.hybrid_affine_character import bounded_affine_positive_roots

    algebra = AffineLieAlgebra(["A", 2, 1])
    roots = bounded_affine_positive_roots(algebra, max_grade=1)

    assert len(roots) == 10
    assert sum(root.is_imaginary for root in roots) == 1
    imaginary_root = next(root for root in roots if root.is_imaginary)
    assert imaginary_root.coordinates == (1, 1, 1)
    assert imaginary_root.multiplicity == 2
    assert all(coordinate >= 0 for root in roots for coordinate in root.coordinates)


@pytest.mark.sage
def test_bounded_affine_verma_partition_low_a1_coefficients():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.hybrid_affine_character import BoundedAffineVermaPartition

    algebra = AffineLieAlgebra(["A", 1, 1])
    partition = BoundedAffineVermaPartition(algebra, max_coordinates=(2, 2))

    assert partition.multiplicity((0, 0)) == 1
    assert partition.multiplicity((0, 1)) == 1
    assert partition.multiplicity((1, 0)) == 1
    assert partition.multiplicity((1, 1)) == 2
    assert partition.multiplicity((-1, 0)) == 0
    assert partition.multiplicity((3, 0)) == 0


@pytest.mark.sage
def test_imaginary_root_partition_uses_finite_rank_multiplicity():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.hybrid_affine_character import BoundedAffineVermaPartition

    algebra = AffineLieAlgebra(["D", 4, 1])
    partition = BoundedAffineVermaPartition(algebra, max_coordinates=(1, 2, 3, 2, 2))

    # The finite D4 partitions contribute 45, while the four independent
    # imaginary-root colors contribute another four copies of delta.
    assert partition.multiplicity((1, 1, 2, 1, 1)) == 49


@pytest.mark.sage
def test_larger_coordinate_box_preserves_lower_partition_coefficients():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.hybrid_affine_character import BoundedAffineVermaPartition

    algebra = AffineLieAlgebra(["A", 2, 1])
    smaller = BoundedAffineVermaPartition(algebra, max_coordinates=(1, 2, 2))
    larger = BoundedAffineVermaPartition(algebra, max_coordinates=(2, 3, 3))

    for coordinates, coefficient in smaller._multiplicities.items():
        assert larger.multiplicity(coordinates) == coefficient


@pytest.mark.sage
def test_partition_rejects_invalid_coordinate_dimensions():
    from sage.all import QQ

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.hybrid_affine_character import BoundedAffineVermaPartition

    algebra = AffineLieAlgebra(["A", 2, 1])
    with pytest.raises(ValueError, match="one entry per affine simple root"):
        BoundedAffineVermaPartition(algebra, max_coordinates=(1, 1))
    with pytest.raises(ValueError, match="bounds must be integers"):
        BoundedAffineVermaPartition(algebra, max_coordinates=(1, 1, QQ(1) / 2))

    partition = BoundedAffineVermaPartition(algebra, max_coordinates=(1, 1, 1))
    with pytest.raises(ValueError, match="one coordinate per affine simple root"):
        partition.multiplicity((0, 0))


@pytest.mark.sage
@pytest.mark.slow
def test_all_d4_order_four_support_partitions_match_recursive_reference():
    from functools import lru_cache

    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.hybrid_affine_character import BoundedAffineVermaPartition

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
    coordinates = tuple(
        candidate.difference_coordinates
        for target in targets
        for candidate in orbit.candidates_above(target)
    )
    bounds = tuple(
        max(query[index] for query in coordinates)
        for index in range(algebra.rank + 1)
    )
    partition = BoundedAffineVermaPartition(algebra, bounds)
    finite_root_lattice = algebra._finite_root_system.root_lattice()
    finite_indices = tuple(finite_root_lattice.index_set())
    finite_positive_roots = tuple(finite_root_lattice.positive_roots())
    marks = tuple(int(algebra.marks[index]) for index in finite_indices)

    def root_coordinates(finite_root, grade):
        coefficients = finite_root.monomial_coefficients()
        return (grade,) + tuple(
            int(coefficients.get(index, 0)) + grade * mark
            for index, mark in zip(finite_indices, marks, strict=True)
        )

    colored_roots = [root_coordinates(root, 0) for root in finite_positive_roots]
    for grade in range(1, bounds[0] + 1):
        colored_roots.extend(
            root_coordinates(root, grade)
            for root in finite_positive_roots + tuple(-root for root in finite_positive_roots)
        )
        colored_roots.extend([(grade,) + tuple(grade * mark for mark in marks)] * algebra.rank)
    colored_roots = tuple(
        root
        for root in colored_roots
        if all(
            root_coordinate <= bound
            for root_coordinate, bound in zip(root, bounds, strict=True)
        )
    )

    @lru_cache(maxsize=None)
    def reference_multiplicity(root_index, remaining):
        if not any(remaining):
            return 1
        if root_index == len(colored_roots):
            return 0
        root = colored_roots[root_index]
        value = reference_multiplicity(root_index + 1, remaining)
        if all(
            root_coordinate <= coordinate
            for root_coordinate, coordinate in zip(root, remaining, strict=True)
        ):
            value += reference_multiplicity(
                root_index,
                tuple(
                    coordinate - root_coordinate
                    for coordinate, root_coordinate in zip(remaining, root, strict=True)
                ),
            )
        return value

    assert len(coordinates) == 238
    assert all(
        partition.multiplicity(query) == reference_multiplicity(0, query)
        for query in coordinates
    )
