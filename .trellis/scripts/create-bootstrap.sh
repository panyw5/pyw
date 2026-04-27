#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common/paths.sh"
source "$SCRIPT_DIR/common/developer.sh"

RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

TASK_NAME="00-bootstrap-guidelines"
PROJECT_TYPE="${1:-research}"

case "$PROJECT_TYPE" in
  research|generic|frontend|backend|fullstack)
    ;;
  *)
    echo -e "${YELLOW}Unknown project type: $PROJECT_TYPE, defaulting to research${NC}"
    PROJECT_TYPE="research"
    ;;
esac

write_prd() {
  local dir="$1"
  cat > "$dir/prd.md" << 'EOF'
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
EOF
}

write_task_json() {
  local dir="$1"
  local developer="$2"
  local today=$(date +%Y-%m-%d)

  cat > "$dir/task.json" << EOF
{
  "id": "$TASK_NAME",
  "name": "Bootstrap Guidelines",
  "description": "Fill project development guidelines for this mathematical physics library",
  "status": "in_progress",
  "dev_type": "docs",
  "priority": "P1",
  "creator": "$developer",
  "assignee": "$developer",
  "createdAt": "$today",
  "completedAt": null,
  "commit": null,
  "subtasks": [
    {"name": "Fill core guidelines", "status": "pending"},
    {"name": "Fill executable test guidelines", "status": "pending"},
    {"name": "Fill mathematical validation guidelines", "status": "pending"},
    {"name": "Add code and paper examples", "status": "pending"}
  ],
  "relatedFiles": [
    ".trellis/spec/core/",
    ".trellis/spec/core-tests/",
    ".trellis/spec/math-tests/",
    ".trellis/spec/guides/"
  ],
  "notes": "First-time setup task adapted for a mathematical physics research library"
}
EOF
}

main() {
  local repo_root=$(get_repo_root)
  local developer=$(get_developer "$repo_root")

  if [[ -z "$developer" ]]; then
    echo -e "${RED}Error: Developer not initialized${NC}"
    echo "Run: ./$DIR_WORKFLOW/$DIR_SCRIPTS/init-developer.sh <your-name>"
    exit 1
  fi

  local tasks_dir=$(get_tasks_dir "$repo_root")
  local task_dir="$tasks_dir/$TASK_NAME"
  local relative_path="$DIR_WORKFLOW/$DIR_TASKS/$TASK_NAME"

  if [[ -d "$task_dir" ]]; then
    echo "$relative_path"
    exit 0
  fi

  mkdir -p "$task_dir"
  write_task_json "$task_dir" "$developer"
  write_prd "$task_dir"
  set_current_task "$relative_path" "$repo_root"

  echo "$relative_path"
}

main "$@"
