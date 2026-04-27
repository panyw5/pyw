"""
调试：检查哪些 Weyl 元素的 Q_tilde 系数为 0
"""

from sage.all import *
import sys

sys.path.insert(0, "/Users/lelouch/pyw")

from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

print("=" * 80)
print("调试：检查哪些 Weyl 元素的 Q_tilde 系数为 0")
print("=" * 80)

alg = AffineLieAlgebra(["D", 4, 1])
kl_character = KazhdanLusztigCharacter(alg)
ωhat = alg.fundamental_weights()
lambda_hat = -2 * ωhat[0]
order = 0

print(f"\nλ̂ = {lambda_hat}")
print(f"order = {order}")

# 准备数据
data = kl_character.prepare_data(lambda_hat, order=order, debug=True)
weyl_list = data.weyl_to_be_summed()

print(f"\nweyl_to_be_summed 项数: {len(weyl_list)}")

# 手动计算每个 Weyl 元素的 Q_tilde 系数
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

zero_count = 0
nonzero_count = 0

print(f"\n计算前 20 个 Weyl 元素的 Q_tilde 系数:")
for i, representative in enumerate(weyl_list[:20]):
    rep_word = tuple(int(j) for j in representative.reduced_word())
    rep_elem = kl_character.kl.weyl_group.from_reduced_word(list(rep_word))
    word_cache[id(rep_elem)] = rep_word

    coefficient = kl_character.kl.affine_bounded_parabolic_Q_tilde_experiment(
        lower_elem,
        rep_elem,
        candidates=all_affine_elems,
        stabilizer_candidates=all_stab_elems,
        at_one=True,
        word_cache=word_cache,
    )

    if coefficient == 0:
        zero_count += 1
        print(f"  {i}: len={representative.length()}, Q̃=0 ❌")
    else:
        nonzero_count += 1
        print(f"  {i}: len={representative.length()}, Q̃={coefficient} ✓")

print(f"\n前 20 项统计:")
print(f"  非零: {nonzero_count}")
print(f"  零: {zero_count}")

print("\n" + "=" * 80)
