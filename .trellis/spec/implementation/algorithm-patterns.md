# FORBIDDEN

- You are **FORBIDDEN** to tell what the user should do or verify: do and verify it yourself, and show it to user
- You are **FORBIDDEN** to answer **BEFORE** you actually run the code user provided or mentioned, and actually see what the error or output is.

## Bounded Affine Verma Partition

The reusable implementation is
`pyw.core.hybrid_affine_character.BoundedAffineVermaPartition`:

```python
partition = BoundedAffineVermaPartition(
    algebra,
    max_coordinates=(max_n_0, ..., max_n_r),
)
multiplicity = partition.multiplicity((n_0, ..., n_r))
```

`max_coordinates` is a componentwise box in the affine simple-root basis. Each
entry must be a nonnegative integer. The table contains the factors for
$\alpha+n\delta$, $-\alpha+n\delta$, and $n\delta$; real roots have
multiplicity one and imaginary roots have finite-rank multiplicity.

Validation matrix:

| Input | Behavior |
|---|---|
| Untwisted affine algebra, integral nonnegative box | Build exact sparse table |
| Query inside box and root cone | Return nonnegative integer |
| Query negative, fractional, or outside box | Return `0` |
| Wrong query dimension | Raise `ValueError` |
| Fractional or negative box | Raise `ValueError` |
| Finite Cartan type | Raise `ValueError` |
| Twisted affine type | Raise `NotImplementedError` |

Required test:

```bash
sage -python -m pytest --no-cov pyw/tests/test_affine_verma_partition.py
```

## Affine Highest-Weight Freudenthal Engine

The bounded engine is implemented in
`pyw.core.hybrid_affine_character`:

```python
tree = build_bounded_weight_tree(lambda_hat, order=order)
multiplicities = compute_bounded_freudenthal_multiplicities(
    tree,
    degenerate_multiplicity=kl_multiplicity,
)
```

Current input contract:

| Input | Behavior |
|---|---|
| Finite-dominant integral affine highest weight | Build bounded domain |
| Nonnegative integral `order` | Include grades `0` through `order` |
| Arbitrary level and highest-weight grade | Use relative grade depth |
| Nonintegral or negative finite Dynkin label | Raise `NotImplementedError` |
| Negative or fractional `order` | Raise `ValueError` |
| Zero Freudenthal denominator | Invoke `degenerate_multiplicity(weight)` |
| Fractional or negative returned multiplicity | Raise `ArithmeticError` |

Only finite-Weyl dominant representatives are recursively computed. The input
contract explicitly requires finite-simple-root integrability. Every
weight above a target is converted to the dominant representative at the same
level and grade before its cached multiplicity is read.

Required test:

```bash
sage -python -m pytest --no-cov pyw/tests/test_hybrid_freudenthal.py
```

## Hybrid Freudenthal-KL Character

The complete bounded character path is
`pyw.core.hybrid_affine_character.KazhdanLusztigFreudenthalCharacter`:

```python
engine = KazhdanLusztigFreudenthalCharacter(algebra)
dominant_multiplicities = engine.multiplicities(lambda_hat, order=order)
explicit_character = engine.character(lambda_hat, order=order)
stats = engine.profile_stats()
```

The zero-denominator branch must apply filters in this order: affine grade,
affine positive-root cone, Verma partition multiplicity, then Q-tilde. One
Q-tilde cache keyed by exact reduced-word tuples is shared by all degenerate
targets in one calculation. `multiplicities()` returns finite-Weyl dominant
weights by relative depth; `character()` expands these multiplicities over the
complete bounded weight tree.

Required test:

```bash
sage -python -m pytest --no-cov pyw/tests/test_hybrid_freudenthal_kl.py
```
