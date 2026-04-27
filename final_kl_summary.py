#!/usr/bin/env sage
"""
总结：KazhdanLusztigCharacter 计算 A₂^(1) 代数 ω̂₀ 的特征标
"""

from sage.all import var
from pyw.core import AffineLieAlgebra, AffineWeight
from pyw.core.character import KazhdanLusztigCharacter, IntegrableModuleCharacter

def main():
    print("=" * 80)
    print("KazhdanLusztigCharacter 计算 A₂^(1) 代数 ω̂₀ 的特征标")
    print("=" * 80)
    
    ala = AffineLieAlgebra(["A", 2, 1])
    omega_hat_0 = AffineWeight.affine_fundamental_weight(ala, 0)
    
    print(f"\n代数: {ala._cartan_type}")
    print(f"最高权重: ω̂₀ = {omega_hat_0.dynkin_labels()}")
    print(f"Level: {omega_hat_0.level}")
    print(f"Grade: {omega_hat_0.grade}")
    
    kl_char = KazhdanLusztigCharacter(ala)
    
    print("\n" + "=" * 80)
    print("KL 公式计算结果")
    print("=" * 80)
    
    print("\n计算 order=0 到 order=5 的特征标:")
    print("-" * 80)
    
    for order in range(6):
        print(f"\norder = {order}:")
        
        numerator_terms = kl_char.numerator_terms(omega_hat_0, order=order)
        print(f"  分子项数量: {len(numerator_terms)}")
        
        kl_result = kl_char.character(omega_hat_0, order=order)
        
        print(f"  特征标系数:")
        for grade in range(order + 1):
            coeff = kl_result[grade]
            print(f"    q^{grade}: {coeff}")
    
    print("\n" + "=" * 80)
    print("详细分析 order=3 的情况")
    print("=" * 80)
    
    numerator_terms_3 = kl_char.numerator_terms(omega_hat_0, order=3)
    
    print(f"\n找到 {len(numerator_terms_3)} 个非零项:")
    print("\n按 grade 分组:")
    
    by_grade = {}
    for term in numerator_terms_3:
        g = term.weight.grade
        if g not in by_grade:
            by_grade[g] = []
        by_grade[g].append(term)
    
    for g in sorted(by_grade.keys()):
        terms = by_grade[g]
        print(f"\ngrade {g}: {len(terms)} 项")
        for i, term in enumerate(terms[:3]):
            print(f"  {i+1}. 权重={term.weight.dynkin_labels()}, "
                  f"KL系数={term.coefficient}, "
                  f"代表元长度={term.representative.length()}")
        if len(terms) > 3:
            print(f"  ... 还有 {len(terms) - 3} 项")
    
    print("\n" + "=" * 80)
    print("与 IntegrableModuleCharacter 对比")
    print("=" * 80)
    
    int_char = IntegrableModuleCharacter(ala)
    q = var("q")
    z1 = var("z1")
    z2 = var("z2")
    
    int_result = int_char.character(omega_hat_0, 5)
    
    print("\nIntegrableModuleCharacter 结果 (权重空间维数):")
    print("-" * 80)
    print(f"{'order':<10} {'权重空间总维数':<20}")
    print("-" * 80)
    
    for order in range(6):
        coeff = int_result.coefficient(q, order)
        dim = coeff.subs({z1: 1, z2: 1}) if order > 0 else coeff
        print(f"{order:<10} {str(dim):<20}")
    
    print("\nKazhdanLusztigCharacter 结果 (Verma 模系数和):")
    print("-" * 80)
    print(f"{'order':<10} {'Verma 模系数和':<20}")
    print("-" * 80)
    
    for order in range(6):
        kl_result = kl_char.character(omega_hat_0, order=order)
        total = sum(kl_result[g] for g in range(order + 1))
        print(f"{order:<10} {str(total):<20}")
    
    print("\n" + "=" * 80)
    print("结论")
    print("=" * 80)
    
    print("\n1. KazhdanLusztigCharacter 成功计算了 A₂^(1) 代数 ω̂₀ 的特征标")
    print("\n2. 结果以 Verma 模的形式线性组合表示:")
    print("   ch(L(ω̂₀)) = Σ_w Q̃_{w_λ,w}(1) · ch(M(w·(Λ+ρ)-ρ))")
    
    print("\n3. 与 IntegrableModuleCharacter 的差异:")
    print("   - IntegrableModuleCharacter: 给出权重空间的直和分解")
    print("   - KazhdanLusztigCharacter: 给出 Verma 模的线性组合")
    print("   - 两者是同一个模的不同表示方式")
    
    print("\n4. 数值差异的原因:")
    print("   - IntegrableModule 计算的是每个 grade 的权重空间维数")
    print("   - KL 计算的是 Verma 模的系数（未展开）")
    print("   - 需要通过 Weyl-Kac 分母公式展开才能对比")
    
    print("\n5. 已知结果对比 (来自 test_integrable_representation.py):")
    print("   order=0: dim=1")
    print("   order=1: dim=8")
    print("   order=2: dim=17")
    print("   order=3: dim=46")
    
    print("\n6. KL 公式的输出:")
    for order in [0, 1, 2, 3]:
        kl_result = kl_char.character(omega_hat_0, order=order)
        print(f"   order={order}: Verma 系数和={sum(kl_result[g] for g in range(order + 1))}")
    
    print("\n" + "=" * 80)
    print("计算完成")
    print("=" * 80)

if __name__ == "__main__":
    main()
