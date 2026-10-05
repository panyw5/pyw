# 迁移自 Trellis 任务 `05-12-D4-level-(-2)`（P1, in_progress, 2026-05-12 起）

## 目标
Consider affine Lie algebra $(\widehat{D}_4)_{-2}$. Consider $\widehat \lambda = -2 \widehat \omega_1$. Compute the character using **Kazhdan-Lusztig method** to $q^0$ order.

## 验收事实
- The character `ch` starts with $1$（q^0 leading term = highest weight state）
- At level $k = -2$ for $D_4$，dual Coxeter number $h^\vee = 6$，so $k + h^\vee = 4$
- $\lambda = -2\omega_1$ corresponds to $\hat\lambda = -2\hat\Lambda_1$（level $= -2$）

## 约束
- use `pyw` to compute. Key classes: `KazhdanLusztigCharacter`, `KazhdanLusztigPolynomials`
- **FORBIDDEN**: 不允许使用 `Algebra.py`（deprecated implementation）
- **FORBIDDEN**: 不允许使用 `IntegrableModuleCharacter`

## 原始资料
- 原 Trellis 目录：`.trellis/tasks/05-12-D4-level-(-2)/`（含 prd.md、run_d4_omega1.py、result_d4_omega1.txt）

## 相关文件
pyw/core/character.py, pyw/core/kazhdan_lusztig.py