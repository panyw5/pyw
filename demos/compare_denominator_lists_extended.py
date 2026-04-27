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
old_den = alg_old.CharacterDen(order)

alg_new = AffineLieAlgebra(cartan_type)
kl = KazhdanLusztigCharacter(alg_new)
new_den = kl.denominator_weight_terms(order)

old_set = set()
for entry in old_den:
    wt = list(entry.keys())[0]
    coeff = int(list(entry.values())[0])
    old_set.add((tuple(wt.to_vector()), coeff))

new_set = set()
for entry in new_den:
    wt = list(entry.keys())[0]
    coeff = int(list(entry.values())[0])
    sage_ext = wt.to_sagemath(extended=True)
    new_set.add((tuple(sage_ext.to_vector()), coeff))

print('old_count', len(old_set))
print('new_count', len(new_set))
print('same', old_set == new_set)
print('only_old', len(old_set - new_set))
print('only_new', len(new_set - old_set))
if old_set != new_set:
    print('sample old-only', sorted(list(old_set - new_set))[:10])
    print('sample new-only', sorted(list(new_set - old_set))[:10])
