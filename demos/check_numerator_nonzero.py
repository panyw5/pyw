"""
检查旧方法 numerator 中非零项的数量
"""

from sage.all import *
import sys
sys.path.insert(0, '/Users/lelouch/pyw/demos')

from MyAlgebra import Alg

print("=" * 80)
print("检查旧方法 numerator 中非零项的数量")
print("=" * 80)

cartanType = ["D", 4, 1]
alg_old = Alg(cartanType, QLoad=True)
order = 0
llambda = - 2 * alg_old.omega[0]

print(f"\nλ = {llambda}")
print(f"order = {order}")

print("\n计算 numerator ...")
numerator_old = alg_old.Kazhdan_Lusztig_numerator(llambda, order)
print(f"numerator 总项数: {len(numerator_old)}")

# 统计非零系数的项
q = var("q")
non_zero_count = 0
zero_count = 0

for term in numerator_old:
    coeff = list(term.values())[0]
    # 代入 q=1
    coeff_at_1 = SR(coeff).subs({q: 1})
    if coeff_at_1 != 0:
        non_zero_count += 1
    else:
        zero_count += 1

print(f"\n非零系数项数: {non_zero_count}")
print(f"零系数项数: {zero_count}")

print(f"\n前 10 项的系数 (代入 q=1):")
for i, term in enumerate(numerator_old[:10]):
    weight = list(term.keys())[0]
    coeff = list(term.values())[0]
    coeff_at_1 = SR(coeff).subs({q: 1})
    print(f"  {i}: weight={weight}, coeff={coeff_at_1}")

print("\n" + "=" * 80)
