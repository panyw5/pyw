# 迁移自 Trellis 任务 `helper-functions`（P2, pending, 2026-04-21 起）

## 目标
Add translator methods for affine weights/roots, list() method for finite Weyl groups, and unify FiniteWeylGroupElement with Sage's native representation.

## 子任务（全部 pending）
- Add translator methods for affine weights and roots（simple root basis, fundamental weight basis）
- Add .list(prefix='s') method to finite_weyl_group for enumerating all elements
- Unify FiniteWeylGroupElement with Sage's native WeylGroup representation

## Notes
Focus on adding convenience methods for basis translation and improving Weyl group API consistency with Sage.

## 相关文件
pyw/core/affine_weight.py, pyw/core/affine_root.py, pyw/core/weyl_group.py, pyw/core/affine_lie_algebra.py, pyw/tests/test_weyl_group.py