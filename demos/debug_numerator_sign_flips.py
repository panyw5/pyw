from sage.all import *
import sys

sys.path.insert(0, "/Users/lelouch/pyw/demos")
sys.path.insert(0, "/Users/lelouch/pyw")
from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

cartan_type = ["D", 4, 1]
order = 0
q = var("q")

alg_old = Alg(cartan_type, QLoad=True)
llambda_old = -2 * alg_old.omega[0]
old_num = alg_old.CharacterNum(llambda_old, order)
old_by_weight = {}
for idx, entry in enumerate(old_num):
    wt = list(entry.keys())[0]
    coeff = SR(list(entry.values())[0]).subs({q: 1})
    old_by_weight[tuple(wt.to_vector())] = (int(coeff), idx, wt)

alg_new = AffineLieAlgebra(cartan_type)
kl_char = KazhdanLusztigCharacter(alg_new)
omega_hat = alg_new.fundamental_weights()
lambda_hat = -2 * omega_hat[0]

data = kl_char.prepare_data(lambda_hat, order=order)
weyl_list = data.weyl_to_be_summed()
num_terms = kl_char.numerator_terms(lambda_hat, order=order)

raw_affine = list(data.W_affine_as_words)
raw_stab = list(data.stabilizer_candidates)
raw_lower = data.w_to_lambda
raw_word_cache = kl_char.kl._build_word_cache([raw_lower] + raw_affine + raw_stab)


def to_elem(w):
    wid = id(w)
    if wid in raw_word_cache:
        return kl_char.kl.weyl_group.from_reduced_word(list(raw_word_cache[wid]))
    return kl_char.kl._ensure_element(w)


all_affine_elems = [to_elem(w) for w in raw_affine]
all_stab_elems = [to_elem(w) for w in raw_stab]
lower_elem = to_elem(raw_lower)

word_cache = {}
for elem, raw in zip(all_affine_elems, raw_affine):
    word_cache[id(elem)] = raw_word_cache[id(raw)]
for elem, raw in zip(all_stab_elems, raw_stab):
    word_cache[id(elem)] = raw_word_cache[id(raw)]
word_cache[id(lower_elem)] = raw_word_cache[id(raw_lower)]

x_max = kl_char.kl._bounded_maximal_representative(lower_elem, stabilizer=all_stab_elems)
word_cache[id(x_max)] = tuple(int(i) for i in x_max.reduced_word())

print("lower_elem =", lower_elem, "word=", word_cache[id(lower_elem)])
print("x_max =", x_max, "word=", word_cache[id(x_max)], "len=", len(word_cache[id(x_max)]))
print("stabilizer =", [tuple(int(i) for i in s.reduced_word()) for s in all_stab_elems])

flip_count = 0
for idx, term in enumerate(num_terms):
    wt = term.weight
    weight_key = tuple(wt.to_sagemath(extended=True).to_vector())
    new_coeff = int(term.coefficient)
    if weight_key not in old_by_weight:
        continue
    old_coeff, old_idx, old_wt = old_by_weight[weight_key]
    if old_coeff == -new_coeff and old_coeff != 0:
        flip_count += 1
        representative = term.representative
        rep_word = tuple(int(i) for i in representative.reduced_word())
        rep_elem = kl_char.kl.weyl_group.from_reduced_word(list(rep_word))
        word_cache[id(rep_elem)] = rep_word
        legacy_val = kl_char.kl.Q_tilde(
            lower_elem,
            rep_elem,
            stabilizer_candidates=all_stab_elems,
            at_one=True,
            word_cache=word_cache,
        )
        bounded_val = kl_char.kl.affine_bounded_parabolic_Q_tilde_experiment(
            lower_elem,
            rep_elem,
            candidates=all_affine_elems,
            stabilizer_candidates=all_stab_elems,
            at_one=True,
            word_cache=word_cache,
        )
        print(
            "\nFLIP",
            flip_count,
            "weight=",
            weight_key,
            "old=",
            old_coeff,
            "new=",
            new_coeff,
            "legacy=",
            legacy_val,
            "bounded=",
            bounded_val,
        )
        print(
            " representative=",
            representative,
            "word=",
            rep_word,
            "len=",
            len(rep_word),
            "old_idx=",
            old_idx,
        )
        coset = []
        for s in all_stab_elems:
            z = rep_elem * s
            z_word = tuple(int(i) for i in z.reduced_word())
            word_cache[id(z)] = z_word
            q_val = kl_char.kl.Q(x_max, z, at_one=True)
            sign = (-1) ** (len(word_cache[id(x_max)]) - len(z_word))
            coset.append((z_word, len(z_word), sign, q_val, sign * q_val))
        print(" coset contributions=", coset)
        if flip_count >= 12:
            break

print("\nflip_count_shown =", flip_count)
