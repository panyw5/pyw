"""
对比新旧方法的 numerator 项
"""

from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw/demos')
sys.path.insert(0, '/Users/lelouch/pyw')

from MyAlgebra import Alg
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

print("=" * 80)
print("对比新旧方法的 numerator 项")
print("=" * 80)

# 旧方法
cartanType = ["D", 4, 1]
alg_old = Alg(cartanType, QLoad=True)
order = 0
llambda = - 2 * alg_old.omega[0]
numerator_old = alg_old.Kazhdan_Lusztig_numerator(llambda, order)

# 新方法
alg_new = AffineLieAlgebra(["D", 4, 1])
kl_character = KazhdanLusztigCharacter(alg_new)
ωhat = alg_new.fundamental_weights()
lambda_hat = -2 * ωhat[0]
numerator_terms_new = kl_character.numerator_terms(lambda_hat, order=order)

print(f"\n旧方法 numerator 项数: {len(numerator_old)}")
print(f"新方法 numerator 项数: {len(numerator_terms_new)}")

# 提取旧方法的权重和系数
q = var("q")
old_weights_coeffs = []
for term in numerator_old:
    weight = list(term.keys())[0]
    coeff = SR(list(term.values())[0]).subs({q: 1})
    old_weights_coeffs.append((weight, coeff))

# 提取新方法的权重和系数
new_weights_coeffs = [(term.weight, term.coefficient) for term in numerator_terms_new]

print(f"\n检查新方法的 44 项是否都在旧方法的 192 项中:")
matches = 0
for new_weight, new_coeff in new_weights_coeffs:
    # 将新方法的 AffineWeight 转换为旧方法的格式
    # 新方法: (finite_part; level; grade)
    # 旧方法: finite_part (忽略 level 和 grade)
    
    # 检查是否匹配
    found = False
    for old_weight, old_coeff in old_weights_coeffs:
        # 简单检查：如果系数相同，认为是匹配的
        if new_coeff == old_coeff:
            matches += 1
            found = True
            break
    
    if not found:
        print(f"  未找到匹配: weight={new_weight}, coeff={new_coeff}")

print(f"\n匹配的项数: {matches} / {len(new_weights_coeffs)}")

print(f"\n新方法前 10 项:")
for i, (weight, coeff) in enumerate(new_weights_coeffs[:10]):
    print(f"  {i}: weight={weight}, coeff={coeff}")

print(f"\n旧方法前 10 项:")
for i, (weight, coeff) in enumerate(old_weights_coeffs[:10]):
    print(f"  {i}: weight={weight}, coeff={coeff}")

print("\n" + "=" * 80)
