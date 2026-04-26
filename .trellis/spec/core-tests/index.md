# Core Tests Guidelines

Use `sage` for testing
```bash
sage -python -m pytest pyw/tests/test_affine_kl_context.py pyw/tests/test_affine_kl_character.py
```

---

## What to Capture

- unit and integration test scope
- interface stability expectations
- regression coverage expectations
- performance or resource baselines when relevant
- what must be tested before code changes are treated as safe
