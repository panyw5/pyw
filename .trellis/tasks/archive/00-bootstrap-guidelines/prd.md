# Bootstrap: Fill Project Development Guidelines

## Purpose

Welcome to Trellis. This bootstrap task exists so the project stops looking like a generic software template and starts reflecting the actual workflow of a mathematical physics Python library.

AI agents use `.trellis/spec/` to understand your conventions. Empty or irrelevant specs make them write plausible but misaligned code.

---

## Your Task

Fill the guideline files based on the real structure and practices of this repository.

### Core Mathematics and Implementation

| File | What to Document |
|------|------------------|
| `.trellis/spec/core/directory-structure.md` | Where core code, demos, refs, docs, and tasks belong |
| `.trellis/spec/core/algebraic-objects.md` | How mathematical objects and public APIs are represented |
| `.trellis/spec/core/algorithm-design.md` | How paper-derived constructions become executable algorithms |
| `.trellis/spec/core/quality-guidelines.md` | Forbidden shortcuts, testing expectations, and review standards |
| `.trellis/spec/core/debugging-and-logging.md` | What to print or log when debugging mathematical code |

### Executable Test Guidelines

| File | What to Document |
|------|------------------|
| `.trellis/spec/core-tests/test-structure.md` | How unit and regression tests are grouped |
| `.trellis/spec/core-tests/regression-and-fixtures.md` | Canonical fixtures and bug-preserving regression cases |
| `.trellis/spec/core-tests/quality-guidelines.md` | What executable correctness tests should prove |

### Mathematical Validation Guidelines

| File | What to Document |
|------|------------------|
| `.trellis/spec/math-tests/reference-validation.md` | How to compare against papers, SageMath, or Wolfram |
| `.trellis/spec/math-tests/invariants-and-identities.md` | Which identities and invariants must be checked |
| `.trellis/spec/math-tests/example-selection.md` | Which examples best reveal mathematical mistakes |

---

## How to Fill Them

1. Look at existing code and tests.
2. Look at existing demos and references.
3. Document the conventions the project already follows.
4. Add examples with real file paths.
5. Record anti-patterns that would create mathematically misleading code.

---

## Completion Checklist

- [ ] Core, executable test, and math validation guidelines filled
- [ ] At least 2-3 real examples in each important guide
- [ ] Anti-patterns documented

When done:

```bash
./.trellis/scripts/task.sh finish
./.trellis/scripts/task.sh archive 00-bootstrap-guidelines
```

---

## Why This Matters

After this task:

1. Trellis will inject research-specific context instead of web-app template assumptions.
2. Future coding sessions will distinguish reusable code, executable tests, and mathematical validation.
3. Future developers and AI agents will onboard faster.
