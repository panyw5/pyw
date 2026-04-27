from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw/demos')
sys.path.insert(0, '/Users/lelouch/pyw')
from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

cartanType = ['D',4,1]
alg_old = Alg(cartanType, QLoad=True)
llambda = -2 * alg_old.omega[0]
old_value = alg_old.Kazhdan_Lusztig_numerator(llambda, 0)
alg_old.Kazhdan_Lusztig_denominator(0)
old_char = alg_old.Kazhdan_Lusztig(0)

alg_new = AffineLieAlgebra(['D',4,1])
kl = KazhdanLusztigCharacter(alg_new)
omega_hat = alg_new.fundamental_weights()
new_char = kl.character(-2*omega_hat[0], order=0)
q = var('q')
coeff0 = new_char.coefficient(q, 0)
b1,b2,b3,b4,q = var('b1 b2 b3 b4 q')
print('old_char =', old_char)
print('new_coeff0 =', coeff0)
print('new_at_ones =', simplify(coeff0.subs({b1:1,b2:1,b3:1,b4:1})))
print('old_minus_new =', simplify(old_char - coeff0.subs({b1:1,b2:1,b3:1,b4:1})))
