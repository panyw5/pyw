# 迁移自 Trellis 任务 `04-25-D6-level-(-4)`（P1, in_progress, 2026-04-25 起）

## 目标
For affine Lie algebra (D6^) at level -4 with lambda_hat = -4 omega_hat_0, compute the character using Kazhdan-Lusztig method and validate low-order constraints.

## 验收事实
- Character starts with 1；q^1 系数必须匹配 D6 adjoint character（66 states）

## 原始资料
- 原 Trellis 目录：`.trellis/tasks/04-25-D6-level-(-4)/`（含 prd.md、stage3~stage6 管线脚本、checkpoint、hybrid_cache/intermediate 产物）

## 进度快照（子任务状态）
- [completed] Persist run context and compute dominant-weight reduction for (D6^)_{-4} vacuum
- [in_progress] Enumerate translations and build affine Weyl-group data in restartable chunks
- [pending] Compute stabilizer, quotient representatives, and weyl_to_sum artifacts with chunked persistence
- [pending] Compute chunked Q_tilde contributions and assemble partial numerator files
- [pending] Compute denominator in partial files and merge final character output
- [pending] Validate q^0 and q^1 constraints, including the 66-state adjoint check

## 相关文件
pyw/core/affine_kl.py, pyw/core/character.py, pyw/core/kazhdan_lusztig.py, pyw/tests/test_affine_kl_character.py, examples/verify_pyw.py