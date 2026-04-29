"""Quick test: Q_tilde with and without at_one=True."""

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
lower = data.w_to_lambda
stabilizer = list(data.stabilizer_candidates)

# Pick one representative
rep = data.quotient_representatives[10]
print(f"Representative: {rep}")
print(f"Stabilizer: {stabilizer}")

# Compare Q_tilde with and without at_one=True
q_tilde_at_one = kl_char.kl.Q_tilde(lower, rep, stabilizer_candidates=stabilizer, at_one=True)
q_tilde_poly = kl_char.kl.Q_tilde(lower, rep, stabilizer_candidates=stabilizer, at_one=False)

print(f"\nQ_tilde(at_one=True):  {q_tilde_at_one}  (type={type(q_tilde_at_one)})")
print(f"Q_tilde(at_one=False): {q_tilde_poly}  (type={type(q_tilde_poly)})")
