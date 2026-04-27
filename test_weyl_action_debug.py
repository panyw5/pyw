#!/usr/bin/env sage
"""
Debug script for Weyl group action issue.

Usage: sage test_weyl_action_debug.py
"""

import sys

sys.path.insert(0, "/Users/lelouch/pyw")

from pyw.core import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight

# Create A2^(1) affine Lie algebra
alg = AffineLieAlgebra(["A", 2, 1])

# Get fundamental weights
Pi = alg.affine_fundamental_weights()
print("=== Affine Fundamental Weights ===")
for i in range(3):
    print(f"Pi[{i}] = {Pi[i]}")
print()

# Get affine Weyl group
W_hat = alg.affine_weyl_group()

# Test simple reflections
print("=== Simple Reflections ===")
for i in range(3):
    si = W_hat.simple_reflection(i)
    print(f"s_{i} = {si}")
    print(f"  finite_part: {si.finite_part}")
    print(f"  translation: {si.translation_vector}")
    print()

# Test action on Pi[0]
print("=== Test Action on Pi[0] ===")
weight = Pi[0]
print(f"weight = {weight}")
print(f"  finite_part: {weight.finite_part}")
print(f"  level: {weight.level}")
print(f"  grade: {weight.grade}")
print()

# Test s_0 action
s0 = W_hat.simple_reflection(0)
print(f"s_0 = {s0}")
result = s0.action(weight)
print(f"s_0.action(Pi[0]) = {result}")
print(f"  finite_part: {result.finite_part}")
print(f"  level: {result.level}")
print(f"  grade: {result.grade}")
print()

# Manual verification
print("=== Manual Verification ===")
print("s_0 should be s_theta * t_{-theta^vee}")
theta = alg.theta_hat().finite_part
print(f"theta = {theta}")
# Get theta_vee from the Weyl group
theta_vee = W_hat.theta_coroot()
print(f"theta^vee = {theta_vee}")
print()

# Step 1: Apply translation t_{-theta^vee}
print("Step 1: Apply translation t_{-theta^vee}")
t_result = alg.translation(-theta_vee, weight)
print(f"t_{{-theta^vee}}(Pi[0]) = {t_result}")
print(f"  finite_part: {t_result.finite_part}")
print(f"  level: {t_result.level}")
print(f"  grade: {t_result.grade}")
print()

# Step 2: Apply s_theta
print("Step 2: Apply s_theta to the translated weight")
# Get the finite Weyl element for s_theta
# theta = alpha_1 + alpha_2 for A2
# s_theta = s_1 * s_2 * s_1 (or s_2 * s_1 * s_2)
W_finite = alg.finite_weyl_group()
s1 = W_finite.simple_reflection(1)
s2 = W_finite.simple_reflection(2)
s_theta = s1 * s2 * s1
print(f"s_theta = {s_theta}")
final_result = s_theta(t_result)
print(f"s_theta(t_{{-theta^vee}}(Pi[0])) = {final_result}")
print(f"  finite_part: {final_result.finite_part}")
print(f"  level: {final_result.level}")
print(f"  grade: {final_result.grade}")
print()

print("=== Comparison ===")
print(f"s_0.action(Pi[0]) = {result}")
print(f"Manual computation = {final_result}")
print(f"Match: {result == final_result}")
