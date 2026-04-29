"""Compare all three: CharacterNum, character_numerator, numerator_terms.

Run with: sage -python demos/compare_all_three.py
"""

import sys

sys.path.insert(0, ".")

from sage.all import *
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight
from pyw.core.character import KazhdanLusztigCharacter


def weight_key(w):
    if hasattr(w, "to_vector"):
        return tuple(w.to_vector())
    if hasattr(w, "to_sagemath"):
        return tuple(w.to_sagemath().to_vector())
    return tuple(w)


def compare_all(algebra_type, lam_expr, order):
    print(f"\n{'=' * 60}")
    print(f"Algebra: {algebra_type},  λ = {lam_expr},  order = {order}")
    print(f"{'=' * 60}")

    ala = AffineLieAlgebra(algebra_type)
    fw = ala.fundamental_weights()
    lam = eval(lam_expr, {"fw": fw, "ala": ala})
    kl_char = KazhdanLusztigCharacter(ala)

    # ── character_numerator (pyw) ──
    print("\n[legacy] Running character_numerator …")
    pyw_result = kl_char.character_numerator(lam, order=order)
    pyw_map = {}
    for entry in pyw_result:
        for weight, coeff in entry.items():
            key = weight_key(weight)
            pyw_map[key] = coeff

    # ── numerator_terms (pyw) ──
    print("\n[numerator_terms] Running numerator_terms …")
    terms_result = kl_char.numerator_terms(lam, order=order)
    terms_map = {}
    for term in terms_result:
        key = weight_key(term.weight)
        terms_map[key] = term.coefficient

    # ── CharacterNum (MyAlgebra.py) ──
    print("\n[MyAlgebra] Running CharacterNum …")
    from MyAlgebra import Alg

    alg = Alg(algebra_type, QLoad=False, WLoad=False)
    llambda = eval(lam_expr, {"fw": alg.omega, "ala": alg})
    legacy_result = alg.CharacterNum(llambda, order=order)
    legacy_map = {}
    for entry in legacy_result:
        for weight, coeff in entry.items():
            key = weight_key(weight)
            legacy_map[key] = coeff

    # ── Compare all three ──
    print(f"\n{'─' * 60}")
    print(
        f"Term counts:  legacy(pyw)={len(pyw_map)},  numerator_terms={len(terms_map)},  CharacterNum={len(legacy_map)}"
    )

    all_keys = set(pyw_map.keys()) | set(terms_map.keys()) | set(legacy_map.keys())

    print(f"\n{'─' * 60}")
    print("Per-weight comparison:")
    print(f"{'weight':<30} {'legacy':>8} {'n_terms':>8} {'CharNum':>8} {'match':>6}")
    print(f"{'─' * 30} {'─' * 8} {'─' * 8} {'─' * 8} {'─' * 6}")

    mismatches_legacy_vs_charnum = 0
    mismatches_nterms_vs_charnum = 0

    for key in sorted(all_keys):
        pyw_c = pyw_map.get(key)
        terms_c = terms_map.get(key)
        legacy_c = legacy_map.get(key)

        pyw_str = str(pyw_c) if pyw_c is not None else "—"
        terms_str = str(terms_c) if terms_c is not None else "—"
        legacy_str = str(legacy_c) if legacy_c is not None else "—"

        # Check if legacy matches CharacterNum
        if pyw_c != legacy_c:
            mismatches_legacy_vs_charnum += 1
        # Check if numerator_terms matches CharacterNum
        if terms_c != legacy_c:
            mismatches_nterms_vs_charnum += 1

        match_str = "✓" if (pyw_c == terms_c == legacy_c) else "✗"
        print(f"{str(key):<30} {pyw_str:>8} {terms_str:>8} {legacy_str:>8} {match_str:>6}")

    print(f"\n{'─' * 60}")
    print(f"legacy vs CharacterNum:     {mismatches_legacy_vs_charnum} mismatch(es)")
    print(f"numerator_terms vs CharacterNum: {mismatches_nterms_vs_charnum} mismatch(es)")

    if mismatches_legacy_vs_charnum == 0:
        print("  ✓ legacy MATCHES CharacterNum")
    if mismatches_nterms_vs_charnum == 0:
        print("  ✓ numerator_terms MATCHES CharacterNum")


if __name__ == "__main__":
    compare_all(["D", 4, 1], "-2 * fw[0]", 2)
