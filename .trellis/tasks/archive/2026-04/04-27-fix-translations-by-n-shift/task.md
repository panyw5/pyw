# 04-27-fix-translations-by-n-shift

## Status
in_progress

## Description
Fix `_translations_by_n_shift_impl` and `_translations_by_n_shift_bnb_impl` functions to only accept one upper bound parameter instead of two (`order` and `max_neg_shift`).

## Problem
The NOTE comment at line 322-327 in `character.py` explains:
- `_translation_neg_shift` calculates which translations make n(weight) change by 0 <= -Δn <= upper_bound
- It should only accept ONE upper bound parameter
- Currently it accepts TWO: `order` and `max_neg_shift`
- The correct upper bound is: `order + n(Λhat + ρhat) - n(λhat)`

## Current Call Sites

### 1. Line 608-613: `_auto_translations`
```python
for translation in _translations_by_n_shift_impl(
    self.algebra,
    affine_weight,
    order=order,
    max_neg_shift=QQ(order),  # Missing n(Λhat+ρhat) - n(λhat)
):
```

### 2. Line 807-812: `_translations_by_n_shift` method
```python
return _translations_by_n_shift_impl(
    self.algebra,
    weight,
    order=order,
    max_neg_shift=max_neg_shift,
)
```

### 3. Line 947-952: `_denominator_candidates`
```python
denominator_translations = _translations_by_n_shift_impl(
    self.algebra,
    rho_hat,
    order=order,
    max_neg_shift=order,  # Missing n(Λhat+ρhat) - n(λhat)
)
```

### 4. Line 1023-1027: `character_weight_list`
```python
translation_elements = self._translations_by_n_shift(
    Lambda_plus_rho,
    order=translation_order,
    max_neg_shift=QQ(translation_order),  # Already correct!
)
```

## Solution
1. Remove `max_neg_shift` parameter from both functions
2. Keep only `order` parameter as the total upper bound
3. Update all callers to pass the correct upper bound

## Acceptance Criteria
- [ ] `_translations_by_n_shift_impl` only has `order` parameter (no `max_neg_shift`)
- [ ] `_translations_by_n_shift_bnb_impl` only has `order` parameter (no `max_neg_shift`)
- [ ] All callers pass the correct upper bound
- [ ] Tests pass
