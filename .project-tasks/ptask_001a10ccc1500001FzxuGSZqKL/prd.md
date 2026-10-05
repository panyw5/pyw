# 迁移自 Trellis 任务 `hybrid-freudenthal-kl-character`（P1, in_progress, 2026-07-16 起，phases 1-6 已全部完成）

## 目标
Use pyw's order-bounded affine Weyl enumeration as a lazy KL orbit index, compute Qtilde only for orbit weights needed at Freudenthal-degenerate targets, and recover the remaining character multiplicities by Freudenthal recursion.

## 进度快照（子任务全部 completed）
- [completed] Extract order-bounded KL orbit preparation without Qtilde evaluation
- [completed] Implement exact affine root coordinates and lazy KL support queries
- [completed] Implement bounded affine roots and Verma partition multiplicities
- [completed] Implement bounded weight tree and Freudenthal recursion
- [completed] Assemble hybrid character and finite-type decomposition
- [completed] Verify D4 level -2 vacuum through q^4 and benchmark cold/warm runs
- [completed] Update mathematical definitions and implementation/test specs

## Notes
- Hybrid path targets complete affine-character coefficients through a requested q order, not the complete explicit KL numerator
- 不专属于 D4 / level -2 / vacuum / 七个固定 targets
- D4 order 4 仅作为 benchmark fixture
- 性能：隔离持久化 Q 与 orbit cache 后，cold D4 order 4 = 50.47s，new-process warm = 5.49s，2112 cached orbit representatives，162/162 Q cache hits
- 最近验证：38 passed

## 相关文件
pyw/core/hybrid_affine_character.py, pyw/core/character.py, pyw/core/kazhdan_lusztig.py, pyw/core/affine_weight.py, pyw/core/affine_lie_algebra.py, pyw/core/weyl_group.py, scripts/negative_level.py, pyw/tests/test_affine_kl_context.py, pyw/tests/test_affine_kl_character.py, pyw/tests/test_hybrid_affine_character.py, pyw/tests/test_affine_verma_partition.py, pyw/tests/test_hybrid_freudenthal.py, pyw/tests/test_hybrid_affine_generality.py