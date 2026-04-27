from sage.all import *
import sys

sys.path.insert(0, "/Users/lelouch/pyw")
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

alg = AffineLieAlgebra(["D", 4, 1])
kl_char = KazhdanLusztigCharacter(alg)
omega_hat = alg.fundamental_weights()
lambda_hat = -2 * omega_hat[0]
order = 0

# 先跑一次 numerator_terms，得到 term.coefficient
terms = kl_char.numerator_terms(lambda_hat, order=order)

# 再独立 prepare_data，避免复用上一步中的局部对象
# 然后每次重算前都清空 KL 相关缓存并不用 word_cache
data = kl_char.prepare_data(lambda_hat, order=order)
raw_stab = list(data.stabilizer_candidates)
raw_lower = data.w_to_lambda
raw_word_cache = kl_char.kl._build_word_cache([raw_lower] + raw_stab)


def to_elem(w):
    wid = id(w)
    if wid in raw_word_cache:
        return kl_char.kl.weyl_group.from_reduced_word(list(raw_word_cache[wid]))
    return kl_char.kl._ensure_element(w)


all_stab_elems = [to_elem(w) for w in raw_stab]
lower_elem = to_elem(raw_lower)

mismatch = []
for i, term in enumerate(terms):
    # 尽量隔离缓存影响
    kl_char.kl._Q_cache.clear()
    kl_char.kl._Q_at_one_cache.clear()
    kl_char.kl._invpol_cache.clear()
    kl_char.kl._P_cache.clear()
    kl_char.kl._legacy_Q_cache.clear()

    rep_word = tuple(int(j) for j in term.representative.reduced_word())
    rep_elem = kl_char.kl.weyl_group.from_reduced_word(list(rep_word))
    recomputed = kl_char.kl.Q_tilde(
        lower_elem,
        rep_elem,
        stabilizer_candidates=all_stab_elems,
        at_one=True,
        word_cache=None,
    )
    if recomputed != term.coefficient:
        mismatch.append(
            (
                i,
                rep_word,
                int(term.coefficient),
                recomputed,
                tuple(term.weight.to_sagemath(extended=True).to_vector()),
            )
        )

print("terms_count", len(terms))
print("mismatch_count_no_cache", len(mismatch))
print("sample_mismatch", mismatch[:20])
