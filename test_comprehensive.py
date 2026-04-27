#!/usr/bin/env sage
"""
Test various affine Weyl group elements to find the bug.
"""

import sys

sys.path.insert(0, "/Users/lelouch/pyw")

from pyw.core import AffineLieAlgebra

# Create A2^(1) affine Lie algebra
alg = AffineLieAlgebra(["A", 2, 1])

# Get fundamental weights
Pi = alg.affine_fundamental_weights()
weight = Pi[0]
print(f"Testing with weight = {weight}")
print()

# Get affine Weyl group
W_hat = alg.affine_weyl_group()

# Generate various affine Weyl group elements
s0 = W_hat.simple_reflection(0)
s1 = W_hat.simple_reflection(1)
s2 = W_hat.simple_reflection(2)

test_cases = [
    ("identity", W_hat.identity()),
    ("s_0", s0),
    ("s_1", s1),
    ("s_2", s2),
    ("s_0 * s_1", s0 * s1),
    ("s_0 * s_2", s0 * s2),
    ("s_1 * s_0", s1 * s0),
    ("s_1 * s_2", s1 * s2),
    ("s_2 * s_0", s2 * s0),
    ("s_2 * s_1", s2 * s1),
    ("s_1 * s_2 * s_1", s1 * s2 * s1),
    ("s_0 * s_1 * s_2", s0 * s1 * s2),
    ("s_0 * s_0", s0 * s0),
    ("s_1 * s_1", s1 * s1),
    ("s_0 * s_1 * s_0", s0 * s1 * s0),
]

print("=== Testing various Weyl group elements ===")
for i, (name, w) in enumerate(test_cases):
    result = w.action(weight)
    print(f"{i:2d}. {name:20s} -> {result}")

    # Check for the bug: (0; 0; 3)
    if result.finite_part == 0 and result.level == 0 and result.grade == 3:
        print(f"    *** BUG FOUND! Got (0; 0; 3) ***")
        print(f"    w = {w}")
        print(f"    finite_part: {w.finite_part}")
        print(f"    translation: {w.translation_vector}")
        break

    # Check for any level change (should never happen!)
    if result.level != weight.level:
        print(f"    *** ERROR: Level changed from {weight.level} to {result.level} ***")
        print(f"    w = {w}")
        print(f"    finite_part: {w.finite_part}")
        print(f"    translation: {w.translation_vector}")
        break

print()
print("=== Testing with different weights ===")
for i in range(3):
    w_test = Pi[i]
    print(f"\nTesting with Pi[{i}] = {w_test}")
    for name, w in [("s_0", s0), ("s_1", s1), ("s_2", s2)]:
        result = w.action(w_test)
        print(f"  {name}(Pi[{i}]) = {result}")
        if result.level != w_test.level:
            print(f"    *** ERROR: Level changed! ***")
