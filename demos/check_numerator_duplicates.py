from sage.all import *
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '/Users/lelouch/pyw/demos')
sys.path.insert(0, '/Users/lelouch/pyw')
from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

cartan_type = ['D',4,1]
order = 0
q = var('q')

alg_old = Alg(cartan_type, QLoad=True)
llambda = -2 * alg_old.omega[0]
old_num = alg_old.CharacterNum(llambda, order)
old_keys = []
old_coeffs_by_weight = defaultdict(list)
for i, entry in enumerate(old_num):
    wt = list(entry.keys())[0]
    coeff = SR(list(entry.values())[0]).subs({q:1})
    key = tuple(wt.to_vector())
    old_keys.append(key)
    old_coeffs_by_weight[key].append((i, int(coeff)))

alg_new = AffineLieAlgebra(cartan_type)
kl = KazhdanLusztigCharacter(alg_new)
omega_hat = alg_new.fundamental_weights()
new_terms = kl.numerator_terms(-2 * omega_hat[0], order=order)
new_keys = []
new_coeffs_by_weight = defaultdict(list)
for i, term in enumerate(new_terms):
    key = tuple(term.weight.to_sagemath(extended=True).to_vector())
    new_keys.append(key)
    new_coeffs_by_weight[key].append((i, int(term.coefficient), tuple(int(j) for j in term.representative.reduced_word())))

old_counter = Counter(old_keys)
new_counter = Counter(new_keys)
old_dups = {k:v for k,v in old_counter.items() if v>1}
new_dups = {k:v for k,v in new_counter.items() if v>1}
print('old_total', len(old_keys), 'old_unique', len(old_counter), 'old_dup_weights', len(old_dups))
print('new_total', len(new_keys), 'new_unique', len(new_counter), 'new_dup_weights', len(new_dups))
if old_dups:
    sample = list(old_dups.items())[:10]
    print('old_dup_sample', sample)
    for key,_ in sample[:5]:
        print(' old_dup_detail', key, old_coeffs_by_weight[key])
if new_dups:
    sample = list(new_dups.items())[:10]
    print('new_dup_sample', sample)
    for key,_ in sample[:5]:
        print(' new_dup_detail', key, new_coeffs_by_weight[key])

shared = set(old_counter) & set(new_counter)
mult_mismatch = [(k, old_counter[k], new_counter[k]) for k in shared if old_counter[k] != new_counter[k]]
print('shared_weights', len(shared))
print('multiplicity_mismatches', len(mult_mismatch))
print('sample_mult_mismatch', mult_mismatch[:10])
