# 迁移自 Trellis 任务 `04-29-port-character-num-to-pyw`（P1, pending, 2026-04-29 起）

## 目标
Analyze demos/MyAlgebra.py CharacterNum function step-by-step and port its numerator business logic as a new function in pyw/core/character.py KazhdanLusztigCharacter class.

## 进度快照（子任务状态）
- [completed] Analyze CharacterNum step-by-step business logic
- [pending] Port CharacterNum logic to new function in KazhdanLusztigCharacter
- [pending] Verify ported function produces identical results to existing numerator_terms

## 背景说明
CharacterNum in MyAlgebra.py is the reference implementation for KL numerator. The existing numerator_terms in KazhdanLusztigCharacter is already its pyw port, but uses prepare_data() delegation. This task creates a new function that mirrors CharacterNum's step ordering exactly.

## 相关文件
demos/MyAlgebra.py, pyw/core/character.py, pyw/core/kazhdan_lusztig.py