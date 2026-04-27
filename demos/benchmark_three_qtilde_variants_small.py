from sage.all import *
import sys
import time

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


def build_case(cartan_type, highest_weight_index, coeff, order):
    alg = AffineLieAlgebra(cartan_type)
    kl_char = KazhdanLusztigCharacter(alg)
    omega_hat = alg.fundamental_weights()
    lambda_hat = coeff * omega_hat[highest_weight_index]
    data = kl_char.prepare_data(lambda_hat, order=order)
    lower = data.w_to_lambda
    stabilizer = list(data.stabilizer_candidates)
    reps = data.weyl_to_be_summed()
    return alg, kl_char, lower, stabilizer, reps


def bench_one(name, runner, loops):
    start = time.perf_counter()
    checksum = 0
    for _ in range(loops):
        vals = runner()
        checksum += sum(int(v) for v in vals)
    elapsed = time.perf_counter() - start
    print(f"{name}: {elapsed:.6f}s checksum={checksum}")
    return elapsed


def run_case(title, cartan_type, highest_weight_index, coeff, order, sample_count, loops):
    alg, kl_char, lower, stabilizer, reps = build_case(
        cartan_type, highest_weight_index, coeff, order
    )
    reps = reps[:sample_count]
    print("\n" + "=" * 100)
    print(title)
    print("sample reps =", len(reps), "stabilizer size =", len(stabilizer), "loops =", loops)

    def run_old():
        return [old_direct(kl_char.kl, alg, lower, rep, stabilizer) for rep in reps]

    def run_legacy():
        return [
            kl_char.kl.Q_tilde(
                lower, rep, stabilizer_candidates=stabilizer, at_one=True, word_cache=None
            )
            for rep in reps
        ]

    def run_bounded():
        return [
            kl_char.kl.affine_bounded_parabolic_Q_tilde_experiment(
                lower,
                rep,
                candidates=reps,
                stabilizer_candidates=stabilizer,
                at_one=True,
                word_cache=None,
            )
            for rep in reps
        ]

    t_old = bench_one("old_direct", run_old, loops)
    t_legacy = bench_one("Q_tilde", run_legacy, loops)
    t_bounded = bench_one("affine_bounded_parabolic_Q_tilde_experiment", run_bounded, loops)
    print(
        "ratios: legacy/old = %.3f, bounded/old = %.3f, bounded/legacy = %.3f"
        % (t_legacy / t_old, t_bounded / t_old, t_bounded / t_legacy)
    )


run_case("A2 fund order3 (all reps)", ["A", 2, 1], 0, 1, 3, sample_count=24, loops=10)
run_case("D4 vacuum order0 (first 24 reps)", ["D", 4, 1], 0, -2, 0, sample_count=24, loops=3)
