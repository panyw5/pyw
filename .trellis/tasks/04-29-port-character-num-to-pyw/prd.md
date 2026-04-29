# background

`demos/MyAlgebra.py` 中的 `CharacterNum` 方法（第 956–1068 行）是 Kazhdan-Lusztig character 公式 numerator 的参考实现。其核心逻辑为：

$$
\text{numerator} = \sum_{w' \in W,\; w_{T\lambda}^{-1} \le w'} \tilde Q(w_{T\lambda}^{-1}, w', W_{\Lambda,0}) \cdot \big[w' \cdot (\Lambda+\rho) - \rho\big]
$$

`pyw/core/character.py` 中已有 `numerator_terms()` 方法，它是 CharacterNum 的 pyw 移植版，但**数据准备通过 `prepare_data()` 委托**，逻辑步骤的编排与 CharacterNum 有差异。

# task: port CharacterNum to pyw

在 `KazhdanLusztigCharacter` 类下新增一个函数，**严格镜像** `CharacterNum` 的步骤顺序，不做任何业务逻辑修改：

## 步骤对照表

| # | CharacterNum (MyAlgebra.py) | 新函数 (pyw) |
|---|---|---|
| 1 | 保存 llambda，确定 Λ 和 wTollambda | 调用 `_find_dominant_Lambda` 或接受外部传入 |
| 2 | 计算 translations：order_base = n(Λ+ρ) - n(λ)，枚举 T | 调用 `_translations_by_n_shift` |
| 3 | 构造 W = W_fin × T | 调用 `_build_W_affine_as_words_direct` |
| 4 | 计算稳定子群 WΛ₀ | 调用 `_legacy_stabilizer_and_quotient_representatives` 或内联 |
| 5 | 计算 dot-轨道 w.(Λ+ρ)-ρ | 循环 W_affine_as_words |
| 6 | 陪集去重，保留最短代表 | 字典去重 |
| 7 | Bruhat 筛选：wTollambda ≤ w' | `lower.bruhat_le(wp)` |
| 8 | 计算 Q̃ 系数，组装 numerator | 调用 `self.kl.Q_tilde(...)` |

## 禁止事项

- **禁止**修改任何业务逻辑（不等式方向、dot-作用公式、Q̃ 参数顺序）
- **禁止**合并或跳过 CharacterNum 中的任何步骤
- **禁止**重命名已有的 `numerator_terms` 方法

## 参考文件

- `demos/MyAlgebra.py` 第 956–1068 行：CharacterNum 原始实现
- `demos/MyAlgebra.py` 第 877–883 行：Qtilde 实现
- `pyw/core/character.py` 第 1163–1242 行：现有 numerator_terms
- `pyw/core/character.py` 第 793–1068 行：KazhdanLusztigCharacter 完整类
