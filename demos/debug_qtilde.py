"""Debug Q_tilde output for a specific case."""

import sys

sys.path.insert(0, ".")

from sage.all import *
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight
from pyw.core.character import KazhdanLusztigCharacter

ala = AffineLieAlgebra(["D", 4, 1])
fw = ala.fundamental_weights()
lam = -2 * fw[0]
kl_char = KazhdanLusztigCharacter(ala)

# Run character_numerator to get the terms
result = kl_char.character_numerator(lam, order=2)

# Check the first few terms
print(f"Total terms: {len(result)}")
print("\nFirst 5 terms:")
for i, entry in enumerate(result[:5]):
    for weight, coeff in entry.items():
        print(f"  [{i}] weight={weight}, coeff={coeff}, type={type(coeff)}")

# Check if any terms have polynomial coefficients
q = var("q")
poly_count = 0
for entry in result:
    for weight, coeff in entry.items():
        if hasattr(coeff, "degree") and coeff.degree(q) > 0:
            poly_count += 1
            if poly_count <= 3:
                print(f"\n  Polynomial coeff: {coeff}")

print(f"\nTerms with polynomial coefficients: {poly_count}")
