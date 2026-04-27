# Pre-commit Hook Setup

This project uses pre-commit hooks to ensure code quality and type annotation correctness.

## Installation

```bash
# Install pre-commit
pip install pre-commit

# Install the git hooks
pre-commit install
```

## What Gets Checked

1. **Code Formatting** (black): Automatically formats Python code
2. **Linting** (ruff): Checks for common Python issues
3. **Type Checking** (mypy): Validates type annotations
4. **TYPE_CHECKING Imports**: Ensures all string type annotations have corresponding imports
5. **General Checks**: Trailing whitespace, file endings, YAML/TOML syntax, etc.

## Manual Run

```bash
# Run on all files
pre-commit run --all-files

# Run on specific files
pre-commit run --files pyw/core/affine_lie_algebra.py

# Run only TYPE_CHECKING check
python scripts/check_type_checking.py pyw/core/*.py
```

## TYPE_CHECKING Import Checker

The custom `check_type_checking.py` script ensures that:
- All string type annotations (e.g., `-> "ClassName"`) have corresponding imports
- These imports are in a `TYPE_CHECKING` block
- This prevents IDE autocomplete issues

Example output:
```
❌ pyw/core/affine_lie_algebra.py
   Missing TYPE_CHECKING imports: ExtendedAffineWeylGroup
   Add to TYPE_CHECKING block:
   
   if TYPE_CHECKING:
       from .module import ExtendedAffineWeylGroup
```

## Skipping Hooks

If you need to skip hooks temporarily:
```bash
git commit --no-verify
```

## Configuration

- Pre-commit config: `.pre-commit-config.yaml`
- Type checking rules: `.github/instructions/python-typing.md`
- Script: `scripts/check_type_checking.py`
