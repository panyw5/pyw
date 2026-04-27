"""
对比旧方法和新方法对同一个 Weyl 元素的 Q_tilde 计算
"""

from sage.all import *
import sys

sys.path.insert(0, "/Users/lelouch/pyw/demos")
sys.path.insert(0, "/Users/lelouch/pyw")

from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

print("=" * 80)
print("对比旧方法和新方法的 Q_tilde 计算")
print("=" * 80)

cartanType = ["D", 4, 1]
alg_old = Alg(cartanType, QLoad=True)
order = 0
llambda = -2 * alg_old.omega[0]
numerator_old = alg_old.Kazhdan_Lusztig_numerator(llambda, order)

alg_new = AffineLieAlgebra(["D", 4, 1])
kl_character = KazhdanLusztigCharacter(alg_new)
ωhat = alg_new.fundamental_weights()
lambda_hat = -2 * ωhat[0]

data = kl_character.prepare_data(lambda_hat, order=order)
weyl_list = data.weyl_to_be_summed()

raw_affine = list(data.W_affine_as_words)
raw_stab = list(data.stabilizer_candidates)
raw_lower = data.w_to_lambda

raw_word_cache = kl_character.kl._build_word_cache([raw_lower] + raw_affine + raw_stab)


def to_elem(w):
    wid = id(w)
    if wid in raw_word_cache:
        return kl_character.kl.weyl_group.from_reduced_word(list(raw_word_cache[wid]))
    return kl_character.kl._ensure_element(w)


all_affine_elems = [to_elem(w) for w in raw_affine]
all_stab_elems = [to_elem(w) for w in raw_stab]
lower_elem = to_elem(raw_lower)

word_cache = {}
for elem, raw in zip(all_affine_elems, raw_affine):
    word_cache[id(elem)] = raw_word_cache[id(raw)]
for elem, raw in zip(all_stab_elems, raw_stab):
    word_cache[id(elem)] = raw_word_cache[id(raw)]
word_cache[id(lower_elem)] = raw_word_cache[id(raw_lower)]

print(f"\n对比前 10 个 Weyl 元素的 Q_tilde 系数:")
print(f"{'Index':<6} {'Length':<8} {'Old Q̃':<15} {'New Q̃':<15} {'Match':<8}")
print("-" * 60)

q = var("q")
for i, representative in enumerate(weyl_list[:10]):
    rep_word = tuple(int(j) for j in representative.reduced_word())
    rep_elem = kl_character.kl.weyl_group.from_reduced_word(list(rep_word))
    word_cache[id(rep_elem)] = rep_word

    new_coeff = kl_character.kl.affine_bounded_parabolic_Q_tilde_experiment(
        lower_elem,
        rep_elem,
        candidates=all_affine_elems,
        stabilizer_candidates=all_stab_elems,
        at_one=True,
        word_cache=word_cache,
    )

    old_coeff = SR(list(numerator_old[i].values())[0]).subs({q: 1})

    match = "✓" if old_coeff == new_coeff else "✗"
    print(
        f"{i:<6} {representative.length():<8} {str(old_coeff):<15} {str(new_coeff):<15} {match:<8}"
    )

print("\n" + "=" * 80)
