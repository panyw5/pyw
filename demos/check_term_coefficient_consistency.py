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

data = kl_char.prepare_data(lambda_hat, order=order)
weyl_list = data.weyl_to_be_summed()
terms = kl_char.numerator_terms(lambda_hat, order=order)

raw_affine = list(data.W_affine_as_words)
raw_stab = list(data.stabilizer_candidates)
raw_lower = data.w_to_lambda
raw_word_cache = kl_char.kl._build_word_cache([raw_lower] + raw_affine + raw_stab)


def to_elem(w):
    wid = id(w)
    if wid in raw_word_cache:
        return kl_char.kl.weyl_group.from_reduced_word(list(raw_word_cache[wid]))
    return kl_char.kl._ensure_element(w)


all_stab_elems = [to_elem(w) for w in raw_stab]
lower_elem = to_elem(raw_lower)
word_cache = {}
for elem, raw in zip(all_stab_elems, raw_stab):
    word_cache[id(elem)] = raw_word_cache[id(raw)]
word_cache[id(lower_elem)] = raw_word_cache[id(raw_lower)]

mismatch = []
for i, term in enumerate(terms):
    rep = term.representative
    rep_word = tuple(int(j) for j in rep.reduced_word())
    rep_elem = kl_char.kl.weyl_group.from_reduced_word(list(rep_word))
    word_cache[id(rep_elem)] = rep_word
    recomputed = kl_char.kl.Q_tilde(
        lower_elem,
        rep_elem,
        stabilizer_candidates=all_stab_elems,
        at_one=True,
        word_cache=word_cache,
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
print("mismatch_count", len(mismatch))
print("sample_mismatch", mismatch[:20])
