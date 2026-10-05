# 迁移自 Trellis 任务 `match-v-flavored-supercharacter-c2-kl`（P2, planning, 2026-08-05 起）

## 目标
Match the V flavored supercharacter to the C2 Kazhdan-Lusztig character（研究型规划任务，原 description 为空，需求需在 brainstorm 中进一步明确）。

## 已有 research 成果（迁移自 `.trellis/tasks/08-05-match-v-flavored-supercharacter-c2-kl/research/c2-kl-flavor.md`）
- 目标代数：`AffineLieAlgebra(["C", 2, 1])`，level -2 vacuum = `-2 * AffineWeight.affine_fundamental_weight(algebra, 0)`
- 现有 eager API：`KazhdanLusztigCharacter(algebra).character(vacuum, order=N)`（pyw/core/character.py）— numerator_q_series / denominator_q_series / taylor prefix
- Fugacity convention：eager KL 实现用 b1..b_rank 单变量，monomial 由 simple-root scalar products 构造
- 关键参考：pyw/tests/test_characters_md_benchmarks.py（flavored q series helper 与 normalization）、pyw/tests/characters.md（character normalization 与 fugacity 定义）、pyw/tests/test_weyl_group_semidirect.py（C2^(1) 构造覆盖）
- 现有 flavored benchmark 覆盖声明在 `.trellis/spec/math-tests/literature-crosscheck.md`

## 原始资料
- 原 Trellis 目录：`.trellis/tasks/08-05-match-v-flavored-supercharacter-c2-kl/`（implement.jsonl / check.jsonl / research/）

## 相关文件
pyw/core/character.py, pyw/core/hybrid_affine_character.py, pyw/tests/test_characters_md_benchmarks.py, pyw/tests/characters.md