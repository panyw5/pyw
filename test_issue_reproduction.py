#!/usr/bin/env sage
"""
Reproduce the issue where w_hat.action(weight) returns (0;0;3).
"""

import sys

sys.path.insert(0, "/Users/lelouch/pyw")

from pyw.core import AffineLieAlgebra

# Create A2^(1) affine Lie algebra
alg = AffineLieAlgebra(["A", 2, 1])

# Get fundamental weights (maybe Pi_u means this?)
Pi_u = alg.affine_fundamental_weights()
print("Pi_u[0] =", Pi_u[0])
print()

# Get affine Weyl group
W_hat = alg.affine_weyl_group()

# Generate some affine Weyl group elements
# Maybe WW_hat is a list of elements?
print("=== Generating affine Weyl group elements ===")

# Method 1: Use simple reflections
elements = []
elements.append(W_hat.identity())  # 0
elements.append(W_hat.simple_reflection(0))  # 1
elements.append(W_hat.simple_reflection(1))  # 2
elements.append(W_hat.simple_reflection(2))  # 3

# Some products
s0 = W_hat.simple_reflection(0)
s1 = W_hat.simple_reflection(1)
s2 = W_hat.simple_reflection(2)

elements.append(s0 * s1)  # 4
elements.append(s0 * s2)  # 5
elements.append(s1 * s0)  # 6
elements.append(s1 * s2)  # 7
elements.append(s2 * s0)  # 8
elements.append(s2 * s1)  # 9
elements.append(s1 * s2 * s1)  # 10

WW_hat = elements

print(f"Generated {len(WW_hat)} elements")
print()

# Test element 10
print("=== Testing WW_hat[10] ===")
w_hat = WW_hat[10]
print(f"w_hat = {w_hat}")
print(f"  finite_part: {w_hat.finite_part}")
print(f"  translation: {w_hat.translation_vector}")
print()

weight = Pi_u[0]
print(f"weight = {weight}")
result = w_hat.action(weight)
print(f"w_hat.action(weight) = {result}")
print()

# Check if this is (0;0;3)
if result.finite_part == 0 and result.level == 0 and result.grade == 3:
    print("ERROR: Got (0;0;3) - this is the bug!")
else:
    print(f"Result: ({result.finite_part}; {result.level}; {result.grade})")
