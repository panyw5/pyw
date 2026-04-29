"""Debug: compare weight deduplication between CharacterNum and pyw."""

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

# Get the data
data = kl_char.prepare_data(lam, order=0)
rho_hat = ala.affine_rho()
rho_sage = rho_hat.to_sagemath(extended=True)
target_sage = (data.Lambda_hat + rho_hat).to_sagemath(extended=True)
W = list(data.W_affine_as_words)

# Compute dot-orbit
print("Computing dot-orbit …")
orbit = [w.action(target_sage) - rho_sage for w in W]

# Check unique weights using different methods
print(f"Total orbit elements: {len(orbit)}")

# Method 1: direct comparison (like CharacterNum)
unique_direct = []
for w in orbit:
    if w not in unique_direct:
        unique_direct.append(w)
print(f"Unique by direct comparison: {len(unique_direct)}")

# Method 2: vector tuple (like pyw)
unique_vector = {}
for w in orbit:
    key = tuple(w.to_vector())
    if key not in unique_vector:
        unique_vector[key] = w
print(f"Unique by vector tuple: {len(unique_vector)}")

# Check if the methods give different results
print(f"\nDifference: {len(unique_vector) - len(unique_direct)}")

# Show some examples of weights that are "different" by vector but "same" by direct comparison
if len(unique_vector) > len(unique_direct):
    print("\nExamples of weights that differ:")
    count = 0
    for i, w1 in enumerate(unique_direct):
        for j, w2 in enumerate(unique_direct):
            if i < j and w1 == w2:
                v1 = tuple(w1.to_vector())
                v2 = tuple(w2.to_vector())
                if v1 != v2:
                    print(f"  w1={w1}, v1={v1}")
                    print(f"  w2={w2}, v2={v2}")
                    count += 1
                    if count >= 3:
                        break
        if count >= 3:
            break
