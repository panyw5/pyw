"""Debug: compare Q output between CharacterNum and pyw."""

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
lower = data.w_to_lambda
stabilizer = list(data.stabilizer_candidates)

# Pick a representative
rep = data.quotient_representatives[10]
print(f"Representative: {rep}")
print(f"lower (w_to_lambda): {lower}")

# Compare Q directly
x_max = kl_char.kl._bounded_maximal_representative(lower, stabilizer=stabilizer)
print(f"x_max: {x_max}")

# pyw Q
q_pyw = kl_char.kl.Q(x_max, rep, at_one=False)
print(f"\npyw Q(x_max, rep, at_one=False): {q_pyw}")
print(f"  type: {type(q_pyw)}")

# CharacterNum Q
from MyAlgebra import Alg

alg = Alg(["D", 4, 1], QLoad=False, WLoad=False)
llambda = -2 * alg.omega[0]
alg.GetLambda(llambda)

# Get the same x_max and rep in CharacterNum's format
# This is tricky because the Weyl group elements might be different objects
# Let me just check what CharacterNum's Q returns for some pairs

print("\nCharacterNum Q examples:")
x = alg.wTollambda
y = alg.W[10]
print(f"  x={x}, y={y}")
q_legacy = alg.Q(x, y)
print(f"  Q(x,y) = {q_legacy}")
