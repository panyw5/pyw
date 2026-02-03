# pyw

**Kac-Moody & W-Algebra Computation Library**

A Python library for computing with Kac-Moody algebras, affine Lie algebras, and W-algebras, with special support for fractional level/weights.

## Features

- **Root Systems**: Finite and affine root systems, (co-)roots
- **Weyl Groups**: Finite, affine, and extended affine Weyl groups
- **Weight Spaces**: Support for fractional weights and levels
- **Admissible Weights**: Kac-Wakimoto admissible weight criteria
- **Principal Admissible Weights**: Construction and verification
- **Embeddings**: $sl_2$-triples, grading, centralizers, Levi subalgebras
- **W-algebras**: Quantum Drinfeld-Sokolov reduction

## Installation

```bash
pip install pyw
```

Or for development:

```bash
git clone https://github.com/panyw5/pyw.git
cd pyw
pip install -e ".[dev,jupyter]"
```

## Requirements

- Python >= 3.10
- SageMath >= 10.0

## Quick Start

```python
from pyw import FractionalLevel, AdmissibleWeight, FractionalWeightSpace

# Create a fractional level k = -h^∨ + p/u
level = FractionalLevel(['A', 2, 1], p=4, u=3)
print(f"Level k = {level.level}")  # k = -3 + 4/3 = -5/3

# Create fractional weights
ws = FractionalWeightSpace(['A', 2, 1])
weight = ws.create_fractional_weight({0: (3, 4), 1: (1, 2)})
print(f"Weight = {weight}")  # 3/4*Lambda[0] + 1/2*Lambda[1]

# Check admissibility
admissible = AdmissibleWeight(['A', 2, 1], weight, level)
print(f"Is admissible: {admissible.is_admissible()}")
```

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
result = trans(cw).action((-4/3) * Lambda[0])
print(f"Result: {result}")
```

## Project Structure

```
pyw/
├── core/          # SageMath wrappers for root systems, Weyl groups, weight spaces
├── fractional/    # Fractional level and admissible weight support
├── embedding/     # Lie algebra embeddings (sl2-triples, Levi subalgebras)
├── walgebra/      # W-algebra constructions
├── algorithms/    # Algorithms from research papers
└── utils/         # LaTeX rendering and visualization
```

## Documentation

See [docs/](docs/) for full documentation.

## References

1. Kac, V. G., Wakimoto, M. "On rationality of W-algebras"
2. Shan, D., Xie, D., Yan, W. "Mirror symmetry for circle compactified 4d N=2 SCFTs"
3. Kac, V. G. "Infinite-dimensional Lie algebras" (3rd ed.)

## License

MIT

## Status

**Version**: 0.1.0 (Alpha)

This project is under active development. API may change.
