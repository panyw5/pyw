"""Compare CharacterNum (MyAlgebra.py) with character_numerator_legacy (pyw).

Run with: sage -python demos/compare_character_num.py
"""

import sys

sys.path.insert(0, ".")

from sage.all import *
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight
from pyw.core.character import KazhdanLusztigCharacter


def weight_key(w):
    """Canonical key for comparing weights across implementations."""
    if hasattr(w, "to_vector"):
        return tuple(w.to_vector())
    if hasattr(w, "to_sagemath"):
        return tuple(w.to_sagemath().to_vector())
    return tuple(w)


def compare_results(algebra_type, lam_expr, order):
    print(f"\n{'=' * 60}")
    print(f"Algebra: {algebra_type},  λ = {lam_expr},  order = {order}")
    print(f"{'=' * 60}")

    # ── pyw ──
    ala = AffineLieAlgebra(algebra_type)
    fw = ala.fundamental_weights()
    lam = eval(lam_expr, {"fw": fw, "ala": ala})
    kl_char = KazhdanLusztigCharacter(ala)

    print("\n[pyw] Running character_numerator_legacy …")
    pyw_result = kl_char.character_numerator_legacy(lam, order=order)
    pyw_map = {}
    for entry in pyw_result:
        for weight, coeff in entry.items():
            key = weight_key(weight)
            pyw_map[key] = coeff

    print(f"[pyw] → {len(pyw_map)} term(s)")
    for key, coeff in sorted(pyw_map.items()):
        print(f"  {key}  →  {coeff}")

    # ── MyAlgebra.py ──
    print("\n[MyAlgebra] Running CharacterNum …")
    try:
        from MyAlgebra import Alg

        alg = Alg(algebra_type, QLoad=False, WLoad=False)
        llambda = eval(lam_expr, {"fw": alg.omega, "ala": alg})
        legacy_result = alg.CharacterNum(llambda, order=order)
        legacy_map = {}
        for entry in legacy_result:
            for weight, coeff in entry.items():
                key = weight_key(weight)
                legacy_map[key] = coeff

        print(f"[MyAlgebra] → {len(legacy_map)} term(s)")
        for key, coeff in sorted(legacy_map.items()):
            print(f"  {key}  →  {coeff}")

        # ── Compare ──
        print(f"\n[COMPARE] pyw terms: {len(pyw_map)},  legacy terms: {len(legacy_map)}")

        all_keys = set(pyw_map.keys()) | set(legacy_map.keys())
        mismatches = 0
        for key in sorted(all_keys):
            pyw_coeff = pyw_map.get(key)
            legacy_coeff = legacy_map.get(key)
            if pyw_coeff != legacy_coeff:
                print(f"  MISMATCH  {key}:  pyw={pyw_coeff},  legacy={legacy_coeff}")
                mismatches += 1

        only_pyw = set(pyw_map.keys()) - set(legacy_map.keys())
        only_legacy = set(legacy_map.keys()) - set(pyw_map.keys())
        if only_pyw:
            print(f"  Only in pyw: {len(only_pyw)} key(s)")
            for k in sorted(only_pyw):
                print(f"    {k}  →  {pyw_map[k]}")
        if only_legacy:
            print(f"  Only in legacy: {len(only_legacy)} key(s)")
            for k in sorted(only_legacy):
                print(f"    {k}  →  {legacy_map[k]}")

        if mismatches == 0 and not only_pyw and not only_legacy:
            print("  ✓ PERFECT MATCH")
        else:
            print(f"  ✗ {mismatches} mismatch(es)")

    except ImportError as e:
        print(f"[MyAlgebra] skipped (import error: {e})")
    except Exception as e:
        print(f"[MyAlgebra] ERROR: {e}")


if __name__ == "__main__":
    # Small test case: A2_1, ω_1, order=1
    compare_results(["A", 2, 1], "fw[1]", 1)
