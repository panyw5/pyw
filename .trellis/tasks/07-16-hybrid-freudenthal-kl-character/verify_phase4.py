from collections import Counter

from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight, affine_weight_key
from pyw.core.hybrid_affine_character import (
    build_bounded_weight_tree,
    classically_dominant_affine_weight,
    compute_bounded_freudenthal_multiplicities,
)

algebra = AffineLieAlgebra(["D", 4, 1])
highest_weight = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
tree = build_bounded_weight_tree(highest_weight, order=4)
degenerate_values = {
    (2, (2, 0, 0, 0)): 2,
    (2, (0, 0, 2, 0)): 2,
    (2, (0, 0, 0, 2)): 2,
    (3, (1, 0, 1, 1)): 8,
    (4, (2, 1, 0, 0)): 16,
    (4, (0, 1, 2, 0)): 16,
    (4, (0, 1, 0, 2)): 16,
}
calls = []


def degenerate_multiplicity(weight):
    calls.append(weight)
    labels = tuple(int(weight.dynkin_labels()[index]) for index in range(1, 5))
    return degenerate_values[(-int(weight.grade), labels)]


multiplicities = compute_bounded_freudenthal_multiplicities(
    tree,
    degenerate_multiplicity=degenerate_multiplicity,
)
dimensions = []
zero_counts = []
for depth in range(5):
    values = [
        multiplicities[
            affine_weight_key(classically_dominant_affine_weight(weight))
        ]
        for weight in tree.weights_at_depth(depth)
    ]
    dimensions.append(sum(values))
    zero_counts.append(values.count(0))

print("counts", [len(tree.weights_at_depth(depth)) for depth in range(5)])
print(
    "dominant",
    Counter(
        -int(tree.weight_by_key[key].grade)
        for key in multiplicities
    ),
)
print(
    "calls",
    len(calls),
    [
        (
            -int(weight.grade),
            tuple(weight.dynkin_labels()[index] for index in range(1, 5)),
        )
        for weight in calls
    ],
)
print("dimensions", dimensions, "zero_counts", zero_counts)
