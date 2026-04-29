"""Direct Q comparison between CharacterNum and pyw."""

import sys

sys.path.insert(0, ".")

from sage.all import *

# ── CharacterNum Q ──
from MyAlgebra import Alg

alg = Alg(["D", 4, 1], QLoad=False, WLoad=False)
llambda = -2 * alg.omega[0]
alg.GetLambda(llambda)

# Get WLambda0
order = 0
order_base = int(alg.Tolambdakn(alg.Lambda + alg.rho)[-1] - alg.Tolambdakn(alg.llambda)[-1])
alg.T = alg.get_translations_by_n_shift(alg.Lambda + alg.rho, order_base + order, order_min=None)
alg.W = alg.GetWeylGroupForqSeries(order=order, T=alg.T)
alg.GetWLambda0(alg.Lambda)

# Get wTollambda
x_legacy = alg.wTollambda
print(f"CharacterNum wTollambda: {x_legacy}")
print(f"CharacterNum WLambda0: {alg.WLambda0}")

# Compute xbar = MaxRep
xbar_legacy = alg.MaxRep(x_legacy, alg.WLambda0)
print(f"CharacterNum xbar: {xbar_legacy}")

# Get a specific Weyl element
# The 10th element in the dot-orbit
LambdaOrbitUnderWeylDot = [(w.action(alg.Lambda + alg.rho) - alg.rho) for w in alg.W]
weightsToBeSummed = []
cosets_legacy = []
for i, weight in enumerate(LambdaOrbitUnderWeylDot):
    if weight not in weightsToBeSummed:
        weightsToBeSummed.append(weight)
        cosets_legacy.append(alg.W[i])
    else:
        ind = weightsToBeSummed.index(weight)
        if alg.W[i].length() < cosets_legacy[ind].length():
            cosets_legacy[ind] = alg.W[i]

WeylToBeSummed_legacy = [wp for wp in cosets_legacy if x_legacy.bruhat_le(wp)]

print(f"\nCharacterNum cosets: {len(cosets_legacy)}")
print(f"CharacterNum WeylToBeSummed: {len(WeylToBeSummed_legacy)}")

# Check Q̃ for first few elements
print("\nCharacterNum Q̃ for first 5 WeylToBeSummed:")
for i, wp in enumerate(WeylToBeSummed_legacy[:5]):
    qt = alg.Qtilde(x_legacy, wp, alg.WLambda0)
    print(f"  [{i}] wp={wp}, Q̃={qt}")

# Check Q directly for specific pairs
print("\nCharacterNum Q for specific pairs:")
for i, wp in enumerate(WeylToBeSummed_legacy[:3]):
    for s in alg.WLambda0:
        z = wp * s
        q_val = alg.Q(xbar_legacy, z)
        print(f"  Q(xbar={xbar_legacy}, z={z}) = {q_val}")
