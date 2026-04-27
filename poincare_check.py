#!/usr/bin/env python3
"""
验证 E_6, m=6 情形的 Poincaré 多项式计算
"""

# 不动流形的计数（从表10）
fixed_manifolds = {
    # (轨道, 维数): 个数
    ("E6(a3)", 3): 1,
    ("E6(a3)", 2): 3,
    ("E6(a3)", 1): 10,
    ("E6(a3)", 0): 16,
    ("D5", 2): 1,
    ("D5", 1): 4,
    ("D5", 0): 15,
    ("E6(a1)", 1): 1,
    ("E6(a1)", 0): 5,
    ("E6", 0): 1,
}

print("=" * 60)
print("E_6 在 m=6 情形的 Poincaré 多项式验证")
print("=" * 60)

# 计算总数
total_by_dim = {}
for (orbit, dim), count in fixed_manifolds.items():
    if dim not in total_by_dim:
        total_by_dim[dim] = 0
    total_by_dim[dim] += count

print("\n按维数统计：")
for dim in sorted(total_by_dim.keys(), reverse=True):
    print(f"  dim={dim}: {total_by_dim[dim]} 个不动流形")

total = sum(total_by_dim.values())
print(f"\n总计: {total} 个不动流形")

# 计算 Poincaré 多项式
# P(t) = Σ (不动流形个数) × t^(2×维数)
print("\nPoincaré 多项式计算：")
print("P(t) = ", end="")
terms = []
for dim in sorted(total_by_dim.keys(), reverse=True):
    count = total_by_dim[dim]
    power = 2 * dim
    if power == 0:
        terms.append(f"{count}")
    else:
        if count == 1:
            terms.append(f"t^{power}")
        else:
            terms.append(f"{count}t^{power}")

print(" + ".join(terms))

# 验证与论文的公式
print("\n论文给出的公式：")
print("P(t) = t^6 + 7t^4 + 27t^2 + 57")

print("\n验证：")
expected = {6: 1, 4: 7, 2: 27, 0: 57}
computed = {2 * dim: count for dim, count in total_by_dim.items()}

match = True
for power in [6, 4, 2, 0]:
    exp_val = expected.get(power, 0)
    comp_val = computed.get(power, 0)
    status = "✓" if exp_val == comp_val else "✗"
    print(f"  t^{power} 系数: 期望={exp_val}, 计算={comp_val} {status}")
    if exp_val != comp_val:
        match = False

if match:
    print("\n✓ 验证通过！计算结果与论文一致。")
else:
    print("\n✗ 验证失败！存在不一致。")

# 特殊表示维数的验证
print("\n" + "=" * 60)
print("特殊表示维数验证")
print("=" * 60)

special_dims = {
    "E6": 1,
    "E6(a1)": 6,
    "D5": 20,
    "E6(a3)": 30,
}

print("\n各轨道的特殊表示维数：")
for orbit, dim in special_dims.items():
    print(f"  Φ_{{{orbit}}}: dim = {dim}")

total_special = sum(special_dims.values())
print(f"\n总计: {total_special}")
print(f"不动流形总数: {total}")

if total_special == total:
    print("\n✓ 猜想4.6验证通过！")
    print("  不动流形个数 = Σ dim(Φ_O)")
else:
    print(f"\n✗ 不一致：{total_special} ≠ {total}")

# 按轨道统计
print("\n" + "=" * 60)
print("按轨道统计不动流形")
print("=" * 60)

orbit_totals = {}
for (orbit, dim), count in fixed_manifolds.items():
    if orbit not in orbit_totals:
        orbit_totals[orbit] = 0
    orbit_totals[orbit] += count

print()
for orbit in ["E6", "E6(a1)", "D5", "E6(a3)"]:
    computed_count = orbit_totals.get(orbit, 0)
    expected_dim = special_dims.get(orbit, 0)
    status = "✓" if computed_count == expected_dim else "✗"
    print(f"  {orbit:12s}: {computed_count:2d} 个不动流形, dim(Φ) = {expected_dim:2d} {status}")
