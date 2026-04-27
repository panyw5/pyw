from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw/demos')
sys.path.insert(0, '/Users/lelouch/pyw')
from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

cartan_type = ['D', 4, 1]
order = 0

alg_old = Alg(cartan_type, QLoad=True)
llambda = -2 * alg_old.omega[0]
old_num = alg_old.CharacterNum(llambda, order)

alg_new = AffineLieAlgebra(cartan_type)
kl = KazhdanLusztigCharacter(alg_new)
omega_hat = alg_new.fundamental_weights()
new_num = kl.numerator_terms(-2 * omega_hat[0], order=order)

old_set = set()
for entry in old_num:
    wt = list(entry.keys())[0]
    coeff = list(entry.values())[0]
    coeff = int(SR(coeff).subs({var('q'): 1})) if hasattr(coeff, 'subs') else int(coeff)
    old_set.add((tuple(wt.to_vector()), coeff))

new_set = set()
for term in new_num:
    wt = term.weight
    coeff = term.coefficient
    coeff = int(SR(coeff).subs({var('q'): 1})) if hasattr(coeff, 'subs') else int(coeff)
    sage_ext = wt.to_sagemath(extended=True)
    new_set.add((tuple(sage_ext.to_vector()), coeff))

print('old_count', len(old_set))
print('new_count', len(new_set))
print('same', old_set == new_set)
print('only_old', len(old_set - new_set))
print('only_new', len(new_set - old_set))
if old_set != new_set:
    print('sample old-only', sorted(list(old_set - new_set))[:20])
    print('sample new-only', sorted(list(new_set - old_set))[:20])
