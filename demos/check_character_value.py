"""
检查新方法的 character 在 q=1, b_i=1 时是否等于 1
"""

from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw')

from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.character import KazhdanLusztigCharacter

print("=" * 80)
print("检查 KazhdanLusztigCharacter.character() 的结果")
print("=" * 80)

alg = AffineLieAlgebra(["D", 4, 1])
kl_character = KazhdanLusztigCharacter(alg)
ωhat = alg.fundamental_weights()
lambda_hat = -2 * ωhat[0]
order = 0

print(f"\nλ̂ = {lambda_hat}")
print(f"order = {order}")

q = var("q")

print("\n计算 character ...")
character_new = kl_character.character(lambda_hat, order=order)

print(f"\ncharacter 类型: {type(character_new)}")
print(f"character 表达式: {character_new}")

# 提取 grade=0 的系数
coeff_0 = character_new.coefficient(q, 0)
print(f"\ngrade=0 的系数:")
print(f"  {coeff_0}")

# 代入 b1=b2=b3=1, q=1
b1, b2, b3, b4, q = var('b1 b2 b3 b4 q')
print(f"\n代入 b1=b2=b3=1:")
try:
    result = coeff_0.subs({b1: 1, b2: 1, b3: 1, b4: 1})
    print(f"  结果 = {result}")
    print(f"  简化后 = {simplify(result)}")
except Exception as e:
    print(f"  错误: {e}")

# 尝试数值计算
print(f"\n数值计算 (b1=b2=b3=1):")
try:
    result_numerical = coeff_0.subs({b1: 1, b2: 1, b3: 1, b4: 1}).n()
    print(f"  结果 = {result_numerical}")
except Exception as e:
    print(f"  错误: {e}")

print("\n" + "=" * 80)
