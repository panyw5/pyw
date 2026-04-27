#!/usr/bin/env python3
"""
演示 Weyl 群作用和椭圆元素的概念
"""

import numpy as np
from itertools import permutations


def apply_permutation(perm, vector):
    """
    对向量应用置换
    perm: 置换，如 (0,2,1) 表示 1→1, 2→3, 3→2
    """
    return np.array([vector[i] for i in perm])


def find_fixed_points(perm, dimension=3):
    """
    找到置换作用下的不动点空间的维数
    对于 A_{n-1} 型，我们在满足 sum=0 的子空间中工作
    """
    n = len(perm)

    identity = np.eye(n)
    perm_matrix = np.array([[1 if perm[i] == j else 0 for j in range(n)] for i in range(n)])

    diff_matrix = perm_matrix - identity

    constraint_vector = np.ones(n)

    augmented = np.vstack([diff_matrix, constraint_vector])

    rank = np.linalg.matrix_rank(augmented)

    fixed_space_dim = n - rank

    return fixed_space_dim


def is_elliptic(perm):
    """
    判断置换是否是椭圆的
    椭圆 <=> 不动点空间只有零向量 <=> 维数为0
    """
    return find_fixed_points(perm) == 0


def cycle_notation(perm):
    """
    将置换转换为循环记号
    """
    n = len(perm)
    visited = [False] * n
    cycles = []

    for i in range(n):
        if not visited[i]:
            cycle = []
            j = i
            while not visited[j]:
                visited[j] = True
                cycle.append(j + 1)
                j = perm[j]
            if len(cycle) > 1:
                cycles.append(tuple(cycle))

    if not cycles:
        return "e"
    return " ".join(f"({' '.join(map(str, c))})" for c in cycles)


print("=" * 70)
print("A_2 型 (SL_3) 的 Weyl 群分析")
print("=" * 70)

print("\nWeyl 群 W = S_3 的所有元素：\n")

perms = list(permutations(range(3)))

elliptic_elements = []

for perm in perms:
    fixed_dim = find_fixed_points(perm)
    is_ell = is_elliptic(perm)
    order = 1
    temp = perm
    while temp != tuple(range(3)):
        temp = tuple(perm[i] for i in temp)
        order += 1

    cycle_str = cycle_notation(perm)
    elliptic_mark = "✓ 椭圆" if is_ell else ""

    print(f"  {cycle_str:12s}  阶数={order}  不动空间维数={fixed_dim}  {elliptic_mark}")

    if is_ell:
        elliptic_elements.append((cycle_str, order))

print(f"\n椭圆元素：{len(elliptic_elements)} 个")
for cycle_str, order in elliptic_elements:
    print(f"  {cycle_str}，阶数 = {order}")

print(f"\n椭圆正则数 (ERN) = {elliptic_elements[0][1] if elliptic_elements else 'N/A'}")

print("\n" + "=" * 70)
print("具体例子：3-循环 (1 2 3) 的作用")
print("=" * 70)

perm_123 = (1, 2, 0)

test_vectors = [
    (1, 0, -1),
    (1, -1, 0),
    (2, -1, -1),
    (1, 1, -2),
]

print("\n测试向量（满足 λ₁ + λ₂ + λ₃ = 0）：\n")

for v in test_vectors:
    v_array = np.array(v)
    w_v = apply_permutation(perm_123, v_array)
    is_fixed = np.allclose(v_array, w_v)

    print(f"  λ = {v}")
    print(f"  w·λ = {tuple(w_v)}")
    print(f"  不动点？ {'是' if is_fixed else '否'}")
    print()

print("结论：(1 2 3) 没有非零不动点，所以是椭圆的。")

print("\n" + "=" * 70)
print("对比：对换 (1 2) 的作用")
print("=" * 70)

perm_12 = (1, 0, 2)

print("\n测试向量：\n")

for v in test_vectors:
    v_array = np.array(v)
    w_v = apply_permutation(perm_12, v_array)
    is_fixed = np.allclose(v_array, w_v)

    print(f"  λ = {v}")
    print(f"  w·λ = {tuple(w_v)}")
    print(f"  不动点？ {'是' if is_fixed else '否'}")
    print()

fixed_example = (1, 1, -2)
print(f"不动点例子：λ = {fixed_example}")
print(f"  (1 2)·λ = {tuple(apply_permutation(perm_12, np.array(fixed_example)))}")
print(f"  λ₁ = λ₂ 时，对换 (1 2) 保持向量不变")
print("\n结论：(1 2) 有非零不动点，所以不是椭圆的。")

print("\n" + "=" * 70)
print("一般规律")
print("=" * 70)

print("""
对于 A_{n-1} 型（SL_n）：

1. Weyl 群 W = S_n（对称群）
2. 作用方式：置换坐标
3. 椭圆元素 = 没有不动点的置换 = n-循环
4. 椭圆正则数 = n（n-循环的阶数）

例如：
- A_2 (SL_3): ERN = 3，椭圆元素是 3-循环
- A_3 (SL_4): ERN = 4，椭圆元素是 4-循环
- A_4 (SL_5): ERN = 5，椭圆元素是 5-循环
""")
