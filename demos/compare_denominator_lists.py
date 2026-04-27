from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw/demos')
sys.path.insert(0, '/Users/lelouch/pyw')
from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

cartan_type = ['D', 4, 1]
order = 0

print('=' * 80)
print('比较 CharacterDen 分母列表')
print('=' * 80)

alg_old = Alg(cartan_type, QLoad=True)
old_den = alg_old.CharacterDen(order)

alg_new = AffineLieAlgebra(cartan_type)
kl = KazhdanLusztigCharacter(alg_new)
new_den = kl.denominator_weight_terms(order)

print(f'旧分母项数: {len(old_den)}')
print(f'新分母项数: {len(new_den)}')

old_pairs = []
for entry in old_den:
    weight = list(entry.keys())[0]
    coeff = list(entry.values())[0]
    old_pairs.append((tuple(weight.to_vector()), int(coeff), weight))
old_pairs.sort(key=lambda item: (item[0], item[1]))

new_pairs = []
for entry in new_den:
    weight = list(entry.keys())[0]
    coeff = list(entry.values())[0]
    sage_weight = weight.to_sagemath(extended=False)
    new_pairs.append((tuple(sage_weight.to_vector()), int(coeff), weight, sage_weight))
new_pairs.sort(key=lambda item: (item[0], item[1]))

old_set = {(vec, coeff) for vec, coeff, _ in old_pairs}
new_set = {(vec, coeff) for vec, coeff, _, _ in new_pairs}

only_old = sorted(old_set - new_set)
only_new = sorted(new_set - old_set)

print(f'完全相同: {old_set == new_set}')
print(f'仅旧版有: {len(only_old)}')
print(f'仅新版有: {len(only_new)}')

if only_old:
    print('\n旧版独有前 10 项:')
    for vec, coeff in only_old[:10]:
        print('  ', vec, coeff)

if only_new:
    print('\n新版独有前 10 项:')
    for vec, coeff in only_new[:10]:
        print('  ', vec, coeff)

print('\n前 10 项逐项对照:')
for i in range(min(10, len(old_pairs), len(new_pairs))):
    old_vec, old_coeff, old_weight = old_pairs[i]
    new_vec, new_coeff, new_weight, new_sage = new_pairs[i]
    print(f'[{i}] old={old_vec},{old_coeff} | new={new_vec},{new_coeff} | same={old_vec == new_vec and old_coeff == new_coeff}')
    print(f'    old weight: {old_weight}')
    print(f'    new weight: {new_weight}')
    print(f'    new sage : {new_sage}')
