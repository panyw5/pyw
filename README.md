# pyw

Distribution name: **pyw-algebra**. The Python import name remains **pyw**.

**Root System & Affine Lie Algebra Computation Library**

A SageMath-based Python library for computing with root systems, affine Lie algebras,
Weyl groups, and fractional levels/weights. This is an experimental research library.
W-algebra constructions are a longer-term goal, not a currently available module.

## Features

- **Root Systems**: Finite and affine root systems, (co-)roots
- **Weyl Groups**: Finite, affine, and extended affine Weyl groups
- **Weight Spaces**: Support for fractional weights and levels
- **Admissible Weights**: Kac-Wakimoto admissible weight criteria
- **Principal Admissible Weights / Nondegeneracy**: Incomplete experimental APIs; not validated for this release
- **Embeddings**: $sl_2$-triples, grading, centralizers, Levi subalgebras
- **Characters**: Experimental integrable, boundary-admissible, and Kazhdan-Lusztig character computations

Support varies by Cartan type and algorithm. An API being present does not imply
that all finite, affine, or twisted types are supported. See the tests and module
docstrings for the currently exercised cases.

Principal-admissible weight construction and nondegeneracy checking are still
under development. Their current tests have known failures; do not rely on
these APIs for verified results. They are retained for experimentation, not
included in the release-preparation validation guarantees.

## Installation

Install SageMath first using the [official installation guide](https://doc.sagemath.org/html/en/installation/).
Use SageMath's Python interpreter, not an unrelated system Python. A PyPI release
is being prepared; until it is published, install from a local checkout:

```bash
git clone https://github.com/panyw5/pyw.git
cd pyw
sage -python -m pip install -e . --no-deps
```

`--no-deps` assumes that a complete SageMath installation already provides the
required `sage` modules. The package metadata declares `sagemath-standard`, but
installing that dependency with ordinary pip can require substantial native
libraries and build tools; it is not a portable one-command SageMath installer.

For development, add the test tools in the same environment:

```bash
sage -python -m pip install pytest pytest-cov
sage -python -m pytest
```

Once a release is available, the corresponding installation command in an
existing SageMath environment will be `sage -python -m pip install pyw-algebra --no-deps`.

## Requirements

- Python >= 3.10
- SageMath >= 10.0

The release-preparation checks currently use SageMath 10.4 with Python 3.12.4
on macOS. Other combinations are not yet verified. Some Kazhdan-Lusztig backends
also require SageMath's optional Coxeter support or a `coxeter` executable; the
executable can be selected with `PYW_COXETER_COMMAND`.

Run Python examples with `sage -python`, or use a SageMath notebook kernel.
In Python files, use `Fraction`, `QQ`, or explicit Sage rationals rather than
floating-point division for exact fractional computations.

## Quick Start

```python
from pyw import FractionalLevel, AdmissibleWeight, FractionalWeightSpace

# Create a fractional level k = -h^∨ + p/u
level = FractionalLevel(['A', 2, 1], p=4, u=3)
print(f"Level k = {level.level}")  # k = -3 + 4/3 = -5/3

# Create a vacuum weight at the same central level
ws = FractionalWeightSpace(['A', 2, 1])
weight = ws.create_fractional_weight({0: level.level})
print(f"Weight = {weight}")  # -5/3*Lambda[0]

# Check admissibility
admissible = AdmissibleWeight(['A', 2, 1], weight, level)
print(f"Is admissible: {admissible.is_admissible()}")
```

`AdmissibleWeight` currently supports finite and untwisted affine types. For
affine input, the weight's central level must match the supplied level. Twisted
affine admissibility is not implemented. To construct other fractional weights,
coefficients may be integers, `Fraction` instances, Sage rationals, or pairs
such as `(3, 4)` for three quarters.

### Weight and Coweight Operations

Weights and coweights live in different mathematical spaces, but you can convert between them:

```python
from pyw.core import AffineLieAlgebra

alg = AffineLieAlgebra(['A', 1, 1])

# Get weights and coweights
Lambda = alg.fundamental_weights()
Lambda_check = alg.fundamental_coweights()

# Convert coweight to weight space for arithmetic
cw = -2 * Lambda_check[1]
cw_as_weight = alg.coweight_to_weight(cw)

# Now can mix with weights
from pyw.core.affine_weight import AffineWeight
w_finite = AffineWeight.affine_fundamental_weight(alg, 0).finite_part
mixed = cw_as_weight + w_finite

# Use in translations with fractional coefficients
trans = alg.extended_affine_weyl_group().translation
from sage.all import QQ
result = trans(cw).action(QQ(-4) / 3 * Lambda[0])
print(f"Result: {result}")
```

## Project Structure

```
pyw/
├── core/          # SageMath wrappers for root systems, Weyl groups, weight spaces
├── fractional/    # Fractional level and admissible weight support
├── embedding/     # Lie algebra embeddings (sl2-triples, Levi subalgebras)
├── algorithms/    # Algorithms from research papers
└── utils/         # Predicates and theta-function helpers
```

## Documentation

See [docs/](docs/) for implementation notes and release guidance. API details
are currently documented in module docstrings; a complete hosted manual is not
yet available.

## References

1. Kac, V. G., Wakimoto, M. "On rationality of W-algebras"
2. Shan, D., Xie, D., Yan, W. "Mirror symmetry for circle compactified 4d N=2 SCFTs"
3. Kac, V. G. "Infinite-dimensional Lie algebras" (3rd ed.)

## License

Licensed under the [MIT License](LICENSE). Copyright (c) 2026 panyw5.

## Status

**Version**: 0.1.0 (Alpha)

This project is under active development. API may change.

## Release checks

See [the release checklist](docs/releasing.md) for build, metadata, and installed-package
verification. Passing a build alone does not establish mathematical correctness.
