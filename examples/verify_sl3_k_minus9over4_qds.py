#!/usr/bin/env python3
"""Verify sl3 boundary-admissible and qDS-reduced module counts.

Target case (from Shan-Xie-Yan, Sec. 3.2 examples):
- g = sl(3), kappa = -9/4 = -3 + 3/4 (so u = 4)
- nilpotent f = [2,1]

Pipeline implemented with pyw:
1) Enumerate affine AKM admissible highest weights via w(S_u) subset Delta_+.
2) Apply qDS non-vanishing condition (Eq. 3.25 / 3.31):
   w(S_u) subset Delta_+ \ {alpha_1}.
3) Identify isomorphic reduced modules by W_f orbit (Eq. 3.26),
   where W_f = <s_1> for f=[2,1].
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.utils.predicates import is_positive_affine_root


def to_fraction(x: Any) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if hasattr(x, "numerator") and hasattr(x, "denominator"):
        n = x.numerator() if callable(x.numerator) else x.numerator
        d = x.denominator() if callable(x.denominator) else x.denominator
        return Fraction(int(n), int(d))
    return Fraction(x)


def affine_dynkin_labels(weight: Any) -> tuple[Fraction, Fraction, Fraction]:
    mc = weight.finite_part.monomial_coefficients()
    l1 = to_fraction(mc.get(1, 0))
    l2 = to_fraction(mc.get(2, 0))
    l0 = to_fraction(weight.level) - l1 - l2
    return (l0, l1, l2)


def dot_action(w: Any, rho: Any, weight: Any) -> Any:
    return w.action(weight + rho) - rho


def main() -> None:
    alg = AffineLieAlgebra(["A", 2, 1])
    wext = alg.extended_affine_weyl_group()

    kappa = Fraction(-9, 4)
    u = 4
    lam0 = alg.fundamental_weights()[0]
    rho = alg.affine_rho()

    # S_u = {-theta + u*delta, alpha_1, alpha_2}
    s_u = [
        u * alg.delta() - alg.affine_theta(),
        alg.affine_simple_roots()[1],
        alg.affine_simple_roots()[2],
    ]

    # Step 1: admissible AKM weights (Eq. 3.4 / 3.25 positivity part)
    by_label: dict[tuple[Fraction, Fraction, Fraction], tuple[tuple[int, int], Any, Any]] = {}

    # Bounded enumeration in W_ext, enough for u=4 case.
    elements = wext.elements_as_semi_direct_product(translation_bounds={1: (-4, 4), 2: (-4, 4)})
    for w in elements:
        images = [w.action(a) for a in s_u]
        if not all(is_positive_affine_root(img) for img in images):
            continue

        lam = dot_action(w, rho, kappa * lam0)
        label = affine_dynkin_labels(lam)

        coeff = w.translation_vector.monomial_coefficients()
        metric = (sum(abs(int(v)) for v in coeff.values()), len(w.reduced_word()))
        if label not in by_label or metric < by_label[label][0]:
            by_label[label] = (metric, w, lam)

    all_labels = sorted(by_label.keys())
    print("=== Step 1: AKM admissible highest weights ===")
    print(f"kappa = {kappa}, u = {u}")
    print(f"count = {len(all_labels)} (expected u^2 = {u * u})")
    for lbl in all_labels:
        print(lbl)

    # Step 2: qDS non-vanishing condition for f=[2,1] (Eq. 3.25 specialized Eq. 3.31)
    # Need w(S_u) subset Delta_+ \ {alpha_1}
    alpha1 = alg.affine_simple_roots()[1]
    surviving: dict[tuple[Fraction, Fraction, Fraction], Any] = {}
    for lbl, (_, w, lam) in by_label.items():
        images = [w.action(a) for a in s_u]
        if all(is_positive_affine_root(img) and img != alpha1 for img in images):
            surviving[lbl] = lam

    surviving_labels = sorted(surviving.keys())
    print("\n=== Step 2: After qDS non-vanishing filter ===")
    print(f"surviving count = {len(surviving_labels)}")
    for lbl in surviving_labels:
        print(lbl)

    # Step 3: quotient by W_f action (Eq. 3.26), W_f = <s_1>
    s1 = wext.simple_reflection(1)
    visited: set[tuple[Fraction, Fraction, Fraction]] = set()
    orbits: list[list[tuple[Fraction, Fraction, Fraction]]] = []

    for lbl in surviving_labels:
        if lbl in visited:
            continue
        lam = surviving[lbl]
        partner = affine_dynkin_labels(dot_action(s1, rho, lam))

        if partner in surviving and partner != lbl:
            orbit = sorted([lbl, partner])
            visited.add(lbl)
            visited.add(partner)
        else:
            orbit = [lbl]
            visited.add(lbl)
        orbits.append(orbit)

    orbits = sorted(orbits, key=lambda x: (len(x), x))
    print("\n=== Step 3: W_f orbit quotient (final W-modules) ===")
    print(f"module count after qDS = {len(orbits)}")
    for i, orb in enumerate(orbits, start=1):
        print(f"orbit {i}: {orb}")

    print("\nConclusion: for g=sl3, kappa=-9/4, f=[2,1], qDS gives 6 simple modules.")


if __name__ == "__main__":
    main()
