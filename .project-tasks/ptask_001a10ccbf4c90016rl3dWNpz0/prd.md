# 迁移自 Trellis 任务 `kazhdan-lusztig-core-pyw`（P1, in_progress, 2026-04-17 起）

## 目标
Build a pyw-native Kazhdan-Lusztig core with stable semantics, bounded affine support, unified cache design, and reproducible validation.

## 原始资料
- 原 Trellis 目录：`.trellis/tasks/kazhdan-lusztig-core-pyw/`（prd.md 含完整需求）

## 进度快照（2026-04-20 notes）
- refs/Kazhdan-Lusztig/Algebra.py 已结构化分解并映射 pyw 复用点
- 语义修正已落地 `pyw/core/kazhdan_lusztig.py`：Q = Weyl-group-element inverse KL，Q_tilde = coset/parabolic inverse KL；affine_bounded_Q 为 canonical，affine_bounded_Q_tilde 为 legacy alias；cache value_kind 更新为 Q_at_one（向后兼容 Q_tilde_at_one）
- 新增 `pyw/core/affine_kl.py`：从 (algebra, lambda_hat, order) 构建 affine KL workflow context
- `KazhdanLusztigCharacter`（pyw/core/character.py）首版已组装 grade-only character
- 44 passed: test_kazhdan_lusztig + test_affine_kl_context + test_affine_kl_character
- 修复 D4^(1) k=-2 legacy 对齐：问题定位在 affine_kl.py 的 extended-weight-lattice 重建而非 Weyl prefix ordering；affine_lie_algebra.py marks 已修正，D4 vacuum 与旧 Algebra.py 对齐至 llambda, rho, llambda+rho 及首批 affine Weyl prefix actions
- AffineLieAlgebra / finite_lie_algebra 增加 dim, num_roots, num_positive_roots, num_negative_roots（45 passed）

## 剩余工作
- 完成 D4 三个 benchmark weights 的 lambda_hat → Lambda_hat regression 收尾
- 将修正后的 GetLambda 等价路径传播到更广的 integral/full-affine workflow

## 相关文件
pyw/core/bruhat.py, pyw/core/kazhdan_lusztig.py, pyw/core/weyl_group.py, pyw/core/affine_kl.py, pyw/core/character.py, pyw/tests/test_kazhdan_lusztig.py, refs/Kazhdan-Lusztig/