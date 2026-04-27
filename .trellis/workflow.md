# Research Workflow

> A Trellis workflow profile for mathematical physics and computational mathematics research.

---

## Quick Start

### Step 0: Initialize Researcher Identity

```bash
python3 ./.trellis/scripts/get_developer.py
python3 ./.trellis/scripts/init_developer.py <your-name>
```

This creates your personal workspace under `.trellis/workspace/<your-name>/`.

### Step 1: Read Research Context

Before making changes, read the profile guidelines in `.trellis/spec/`:

- `definitions/` for mathematical and physical definitions, symbols, and canonical forms
- `implementation/` for data representations, algorithm patterns, module boundaries, and caching rules
- `core-tests/` for programming-project tests, regression expectations, and interface coverage
- `math-tests/` for mathematical benchmark cases, identities, and literature cross-checks
- `phy-tests/` for physical observables, symmetries, limits, and consistency expectations
- `guides/` for thinking guides that connect definitions, implementation, and validation

### Step 2: Start from a Research Task

Create or select a task under `.trellis/tasks/`.
Each task should center on a mathematical question, construction, algorithm, or computation experiment.

---

## Research Task Lifecycle

Use this mental flow when working:

1. **Formulate** the problem and notation.
2. **Construct** the mathematical objects and helper routines.
3. **Experiment** with examples, parameter sweeps, or symbolic calculations.
4. **Verify** against known cases, identities, symmetries, and literature.
5. **Synthesize** the conclusion, failed attempts, and open questions.
6. **Archive** the session so the next pass starts with real context.

---

## Completion Standard

Treat work as complete only when the relevant items are true:

- notation is consistent with the project conventions
- key benchmark cases were checked
- output is reproducible
- failure cases or open questions are recorded
- conclusions distinguish proved results from computational evidence

When a task includes package code, also run the project tests and type checks.

---

## Core Principle

Trellis is here to preserve research context across long-running investigations.
Use specs to make definitions and conventions explicit, tasks to organize lines of attack, and workspace journals to remember what was tried and what remains open.
