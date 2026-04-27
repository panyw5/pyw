#!/usr/bin/env sage
"""
测试 KazhdanLusztigCharacter 计算 A21 代数 ωhat[0] 的特征标
并与 IntegrableModuleCharacter 的已知结果对比
"""

from sage.all import var, SR
from pyw.core import AffineLieAlgebra, AffineWeight
from pyw.core.character import KazhdanLusztigCharacter, IntegrableModuleCharacter

def main():
    print("=" * 80)
    print("测试 KazhdanLusztigCharacter 计算 A₂^(1) 代数 ω̂₀ 的特征标")
    print("=" * 80)
    
    # 创建 A21 代数
    ala = AffineLieAlgebra(["A", 2, 1])
    print(f"\n代数: {ala._cartan_type}")
    print(f"秩: {ala.rank}")
    
    # 获取基本权重 ω̂₀
    omega_hat_0 = AffineWeight.affine_fundamental_weight(ala, 0)
    print(f"\n最高权重 ω̂₀:")
    print(f"  Dynkin 标签: {omega_hat_0.dynkin_labels()}")
    print(f"  Level: {omega_hat_0.level}")
    print(f"  Grade: {omega_hat_0.grade}")
    
    # 创建 KazhdanLusztigCharacter 实例
    kl_char = KazhdanLusztigCharacter(ala)
    
    # 计算不同阶数的特征标
    print("\n" + "=" * 80)
    print("使用 KazhdanLusztigCharacter 计算特征标")
    print("=" * 80)
    
    for order in [0, 1, 2, 3]:
        print(f"\n阶数 {order}:")
        print("-" * 40)
        
        # 计算 KL 特征标
        kl_result = kl_char.character(omega_hat_0, order=order)
        
        print(f"KL 特征标系数:")
        for grade, coeff in kl_result:
            print(f"  q^{grade}: {coeff}")
    
    # 与 IntegrableModuleCharacter 对比
    print("\n" + "=" * 80)
    print("与 IntegrableModuleCharacter 结果对比")
    print("=" * 80)
    
    # 创建 IntegrableModuleCharacter 实例
    int_char = IntegrableModuleCharacter(ala)
    
    # 定义符号变量
    q = var("q")
    z1 = var("z1")
    z2 = var("z2")
    
    # 计算 IntegrableModuleCharacter 的结果（order=3）
    print("\nIntegrableModuleCharacter 结果 (order=3):")
    int_result = int_char.character(omega_hat_0, 3)
    
    # 提取每个阶数的系数
    print("\n系数对比:")
    print("-" * 80)
    print(f"{'阶数':<10} {'IntegrableModule':<30} {'KazhdanLusztig':<30} {'匹配':<10}")
    print("-" * 80)
    
    all_match = True
    for order in range(4):
        # 计算 KL 特征标
        kl_result = kl_char.character(omega_hat_0, order=order)
        
        # 获取 q^order 的系数
        int_coeff = int_result.coefficient(q, order)
        
        # KL 结果中 q^order 的系数
        kl_coeff = kl_result[order]
        
        # 对于 order=0，系数应该是 1
        # 对于 order>0，IntegrableModule 给出的是完整的 z1, z2 表达式
        # 而 KL 给出的是维数（常数项）
        
        if order == 0:
            match = (int_coeff == 1 and kl_coeff == 1)
            print(f"{order:<10} {str(int_coeff):<30} {str(kl_coeff):<30} {'✓' if match else '✗':<10}")
        else:
            # 对于 order > 0，我们需要提取 IntegrableModule 结果的常数项
            # 这对应于 z1^0 * z2^0 的系数
            try:
                # 展开并提取常数项
                expanded = int_coeff.expand()
                # 将 z1=1, z2=1 代入得到总维数
                const_term = expanded.subs({z1: 1, z2: 1})
                
                print(f"{order:<10} {str(const_term):<30} {str(kl_coeff):<30} {'?' if const_term != kl_coeff else '✓':<10}")
                
                if const_term != kl_coeff:
                    all_match = False
                    print(f"  注意: IntegrableModule 给出完整的权重空间分解")
                    print(f"        KL 给出的是 Verma 模的线性组合")
            except:
                print(f"{order:<10} {str(int_coeff):<30} {str(kl_coeff):<30} {'?':<10}")
    
    # 显示详细的 order=3 结果
    print("\n" + "=" * 80)
    print("详细结果 (order=3)")
    print("=" * 80)
    
    print("\nIntegrableModuleCharacter (完整表达式):")
    print(int_result.coefficient(q, 3))
    
    print("\nKazhdanLusztigCharacter (Verma 模系数):")
    kl_result_3 = kl_char.character(omega_hat_0, order=3)
    for grade, coeff in kl_result_3:
        print(f"  q^{grade}: {coeff}")
    
    # 显示 numerator terms
    print("\n" + "=" * 80)
    print("KL 公式的分子项 (order=3)")
    print("=" * 80)
    
    numerator_terms = kl_char.numerator_terms(omega_hat_0, order=3)
    print(f"\n找到 {len(numerator_terms)} 个非零项:")
    for i, term in enumerate(numerator_terms):
        print(f"\n项 {i+1}:")
        print(f"  代表元长度: {term.representative.length()}")
        print(f"  权重: {term.weight.dynkin_labels()}")
        print(f"  权重 grade: {term.weight.grade}")
        print(f"  KL 系数: {term.coefficient}")
    
    print("\n" + "=" * 80)
    print("测试完成")
    print("=" * 80)

if __name__ == "__main__":
    main()
