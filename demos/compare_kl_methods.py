"""
对比 MyAlgebra.py 和 KazhdanLusztigCharacter 的计算结果

测试用例: D4^(1), λ = -2 ω̂₀, order = 0
"""

from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw/demos')
sys.path.insert(0, '/Users/lelouch/pyw')

from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter
from pyw.core.affine_weight import AffineWeight

print("=" * 80)
print("测试用例: D4^(1), λ = -2 ω̂₀, order = 0")
print("=" * 80)

# ============================================================================
# 方法 1: MyAlgebra.py (旧方法)
# ============================================================================
print("\n[方法 1] MyAlgebra.Alg")
print("-" * 80)

cartanType = ["D", 4, 1]
alg_old = Alg(cartanType, QLoad=True)
order = 0
llambda = - 2 * alg_old.omega[0]

print(f"λ = {llambda}")
print(f"order = {order}")

print("\n计算 numerator ...")
numerator_old = alg_old.Kazhdan_Lusztig_numerator(llambda, order)
print(f"numerator 项数: {len(numerator_old)}")
print(f"numerator 前 3 项:")
for i, term in enumerate(numerator_old[:3]):
    print(f"  {i}: {term}")

print("\n计算 denominator ...")
denominator_old = alg_old.Kazhdan_Lusztig_denominator(order)
print(f"denominator 项数: {len(denominator_old)}")
print(f"denominator 前 3 项:")
for i, term in enumerate(denominator_old[:3]):
    print(f"  {i}: {term}")

print("\n计算最终 character ...")
character_old = alg_old.Kazhdan_Lusztig(order)
print(f"character = {character_old}")

# ============================================================================
# 方法 2: KazhdanLusztigCharacter (新方法)
# ============================================================================
print("\n\n[方法 2] KazhdanLusztigCharacter")
print("-" * 80)

alg_new = AffineLieAlgebra(["D", 4, 1])
kl_character = KazhdanLusztigCharacter(alg_new)
ωhat = alg_new.fundamental_weights()
lambda_hat = -2 * ωhat[0]

print(f"λ̂ = {lambda_hat}")
print(f"order = {order}")

print("\n计算 numerator terms ...")
numerator_terms_new = kl_character.numerator_terms(lambda_hat, order=order)
print(f"numerator 项数: {len(numerator_terms_new)}")
print(f"numerator 前 3 项:")
for i, term in enumerate(numerator_terms_new[:3]):
    print(f"  {i}: weight={term.weight}, coeff={term.coefficient}")

print("\n计算 character ...")
character_new = kl_character.character(lambda_hat, order=order)
print(f"character = {character_new}")

# ============================================================================
# 对比结果
# ============================================================================
print("\n\n[对比结果]")
print("=" * 80)

print(f"\n1. numerator 项数对比:")
print(f"   旧方法: {len(numerator_old)}")
print(f"   新方法: {len(numerator_terms_new)}")

print(f"\n2. denominator 项数对比:")
print(f"   旧方法: {len(denominator_old)}")
print(f"   新方法: (新方法在 character() 内部计算)")

print(f"\n3. 最终 character 对比:")
print(f"   旧方法: {character_old}")
print(f"   新方法: {character_new}")

# 检查 numerator 的权重是否相同
print(f"\n4. numerator 权重对比:")
print("   旧方法的权重:")
for i, term in enumerate(numerator_old[:5]):
    weight = list(term.keys())[0]
    print(f"     {i}: {weight}")

print("   新方法的权重:")
for i, term in enumerate(numerator_terms_new[:5]):
    print(f"     {i}: {term.weight}")

print("\n" + "=" * 80)
print("对比完成")
print("=" * 80)
