from sage.all import *
import sys

sys.path.insert(0, "/Users/lelouch/pyw")
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter


def old_direct(kl_obj, algebra, x_min, y_min, subgroup):
    group_sage = algebra.affine_weyl_group_sage()
    x_min_sage = group_sage.from_reduced_word([int(i) for i in x_min.reduced_word()])
    y_min_sage = group_sage.from_reduced_word([int(i) for i in y_min.reduced_word()])
    subgroup_list = [
        group_sage.from_reduced_word([int(i) for i in s.reduced_word()]) for s in subgroup
    ]
    x_max = max([x_min_sage * s for s in subgroup_list], key=lambda w: int(w.length()))
    x_max_elem = kl_obj.weyl_group.from_reduced_word([int(i) for i in x_max.reduced_word()])

    result = 0
    for stabilizer_element in subgroup_list:
        coset_element = y_min_sage * stabilizer_element
        coset_elem = kl_obj.weyl_group.from_reduced_word(
            [int(i) for i in coset_element.reduced_word()]
        )
        q_value = kl_obj.Q(x_max_elem, coset_elem, at_one=True)
        sign = int((-1) ** (int(x_max.length()) - int(coset_element.length())))
        result += sign * q_value
    if hasattr(result, "full_simplify"):
        result = result.full_simplify()
    return result


def build_sample(cartan_type, highest_weight_index, coeff, order):
    alg = AffineLieAlgebra(cartan_type)
    kl_char = KazhdanLusztigCharacter(alg)
    omega_hat = alg.fundamental_weights()
    lambda_hat = coeff * omega_hat[highest_weight_index]
    data = kl_char.prepare_data(lambda_hat, order=order)
    lower = data.w_to_lambda
    stabilizer = list(data.stabilizer_candidates)
    reps = data.weyl_to_be_summed()
    return alg, kl_char, lower, stabilizer, reps


def compare_case(name, cartan_type, highest_weight_index, coeff, order, sample_limit=20):
    alg, kl_char, lower, stabilizer, reps = build_sample(
        cartan_type, highest_weight_index, coeff, order
    )
    total = len(reps)
    print("\n" + "=" * 100)
    print("CASE", name)
    print(
        "cartan_type =",
        cartan_type,
        "lambda coeff/index =",
        coeff,
        highest_weight_index,
        "order =",
        order,
    )
    print("total reps =", total)
    print(
        "stabilizer size =",
        len(stabilizer),
        "stabilizer words =",
        [tuple(int(i) for i in s.reduced_word()) for s in stabilizer],
    )

    mism_old_vs_legacy = []
    mism_old_vs_bounded = []
    mism_legacy_vs_bounded = []

    for idx, rep in enumerate(reps):
        old_val = old_direct(kl_char.kl, alg, lower, rep, stabilizer)
        legacy_val = kl_char.kl.Q_tilde(
            lower, rep, stabilizer_candidates=stabilizer, at_one=True, word_cache=None
        )
        bounded_val = kl_char.kl.affine_bounded_parabolic_Q_tilde_experiment(
            lower,
            rep,
            candidates=reps,
            stabilizer_candidates=stabilizer,
            at_one=True,
            word_cache=None,
        )
        row = (idx, tuple(int(i) for i in rep.reduced_word()), old_val, legacy_val, bounded_val)
        if old_val != legacy_val:
            mism_old_vs_legacy.append(row)
        if old_val != bounded_val:
            mism_old_vs_bounded.append(row)
        if legacy_val != bounded_val:
            mism_legacy_vs_bounded.append(row)

    print("old vs legacy mismatches =", len(mism_old_vs_legacy))
    print("old vs bounded mismatches =", len(mism_old_vs_bounded))
    print("legacy vs bounded mismatches =", len(mism_legacy_vs_bounded))
    if mism_old_vs_legacy:
        print("sample old vs legacy:")
        for row in mism_old_vs_legacy[:sample_limit]:
            print(" ", row)
    if mism_old_vs_bounded:
        print("sample old vs bounded:")
        for row in mism_old_vs_bounded[:sample_limit]:
            print(" ", row)
    if mism_legacy_vs_bounded:
        print("sample legacy vs bounded:")
        for row in mism_legacy_vs_bounded[:sample_limit]:
            print(" ", row)


compare_case("D4 vacuum order0", ["D", 4, 1], 0, -2, 0, sample_limit=12)
compare_case("A2 fund order3", ["A", 2, 1], 0, 1, 3, sample_limit=12)
