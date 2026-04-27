#!/usr/bin/env sage
"""
Kazhdan-Lusztig Character Computation Demo

This script demonstrates how to compute characters of admissible modules
using the Kazhdan-Lusztig formula implemented in pyw.

Example: A₂^(1) at fractional level k = -3 + 4/3 = -5/3

Usage:
    sage kazhdan_lusztig_demo.py

Requirements:
    - SageMath >= 10.0
    - pyw package installed
"""

from pyw.core import AffineLieAlgebra, AffineWeight
from pyw.core.bruhat import BruhatOrder
from pyw.core.character import KazhdanLusztigCharacter
from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials
from pyw.fractional import FractionalLevel

from sage.all import WeylGroup, var


def demo_bruhat_order():
    """Demonstrate Bruhat order computations."""
    print("=" * 60)
    print("Demo 1: Bruhat Order for A₂")
    print("=" * 60)

    W = WeylGroup(["A", 2])
    bruhat = BruhatOrder(W)

    e = W.one()
    s1 = W.simple_reflection(1)
    s2 = W.simple_reflection(2)
    w0 = W.long_element()

    print(f"\nWeyl group: {W}")
    print(f"Identity: {e}")
    print(f"s₁: {s1}")
    print(f"s₂: {s2}")
    print(f"Longest element w₀: {w0}")

    print("\nBruhat order comparisons:")
    print(f"  e ≤ s₁: {bruhat.le(e, s1)}")
    print(f"  s₁ ≤ s₁s₂: {bruhat.le(s1, s1 * s2)}")
    print(f"  s₁ ≤ w₀: {bruhat.le(s1, w0)}")

    print("\nBruhat interval [e, w₀]:")
    interval = bruhat.interval(e, w0)
    for w in interval:
        print(f"  ℓ={bruhat.length(w)}: {w}")


def demo_kl_polynomials():
    """Demonstrate KL polynomial computations."""
    print("\n" + "=" * 60)
    print("Demo 2: Kazhdan-Lusztig Polynomials for A₂")
    print("=" * 60)

    W = WeylGroup(["A", 2])
    kl = KazhdanLusztigPolynomials(W)

    e = W.one()
    s1 = W.simple_reflection(1)
    s2 = W.simple_reflection(2)
    w0 = W.long_element()

    print("\nKL polynomials P_{x,y}(1):")
    print(f"  P_{{e,e}}(1) = {kl.P(e, e, at_one=True)}")
    print(f"  P_{{e,s₁}}(1) = {kl.P(e, s1, at_one=True)}")
    print(f"  P_{{e,w₀}}(1) = {kl.P(e, w0, at_one=True)}")

    print("\nFinite-group experimental inverse KL polynomials Q̃_{x,y}(1):")
    print(f"  Q̃_{{e,e}}(1) = {kl.Q_tilde_experiment(e, e, at_one=True)}")
    print(f"  Q̃_{{e,s₁}}(1) = {kl.Q_tilde_experiment(e, s1, at_one=True)}")


def demo_affine_character():
    """Demonstrate the supported affine character entrypoint."""
    print("\n" + "=" * 60)
    print("Demo 3: Kazhdan-Lusztig Character for A₂^(1)")
    print("=" * 60)

    ala = AffineLieAlgebra(["A", 2, 1])
    kl_char = KazhdanLusztigCharacter(ala)
    q = var("q")
    z1 = var("z1")
    z2 = var("z2")
    lam = AffineWeight.affine_fundamental_weight(ala, 0)

    print(f"\nAlgebra: {ala._cartan_type}")
    print("\nComputing character(Λ₀) up to order 2...")

    result = kl_char.character(lam, order=2)

    print(f"Character: {result}")
    print("\nq-coefficients:")
    for grade in range(3):
        print(f"  q^{grade}: {result.coefficient(q, grade)}")
    print(f"\nSpecialized at z1=z2=1: {result.subs({z1: 1, z2: 1})}")


if __name__ == "__main__":
    print("Kazhdan-Lusztig Character Computation Demo")
    print("=" * 60)

    demo_bruhat_order()
    demo_kl_polynomials()
    demo_affine_character()

    print("\n" + "=" * 60)
    print("Demo complete!")
    print("=" * 60)
