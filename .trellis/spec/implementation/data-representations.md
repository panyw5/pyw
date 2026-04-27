# Naming convention

Naming of variables, methods and properties should be self-explanatory, easy to understand for a human researcher or student. Key rules to follow:
- **DIRECT**: explicitly state "what it is"
- **AVOID_ABBRIVIATION**: unless it is well-known, or human name initials, avoid abbreviations.
- **DESCRIPTIVE SUFFIX**: When several variables refer to the same mathematical object but in different representations, or different origin, add **descriptive suffixes**
- **CONSISTENCY**: Naming of variables pointing to the same concept should be **consistent** across the codebase
- **STANDARD MATH**: Use standard mathematical notation: `W_hat` or `W_aff` for $\widehat W$, `omega_hat` for $\widehat \omega$

- Use **PascalCase** for class names, and **lowercase_with_underscores** (snake_case) for variables, methods and properties.


## DO and DONT examples:

```python
rs = alg._finite_root_system()                  # DON'T
finite_root_system = alg._finite_root_system()  # DO
```

```python

semidirect = alg.affine_weyl_group()          # DON'T
affine_weyl_group = alg.affine_weyl_group()   # DO
affine_weyl_group_as_semidirect = alg.affine_weyl_group()   # DO
```

```python
class KazhdanLusztigPolynomials:
  ...

  # DON'T: what these elements belong to is unclear, origin unclear
  def affine_bounded_elements(self, ...):
    r"""Enumerate a bounded affine subset and coerce it into ``self.weyl_group``.
    """

  # DO: elements are clearly from affine weyl group, and from bounded translations
  def affine_weyl_elements_from_bounded_translations(self, ...):
    r"""Enumerate a bounded affine subset and coerce it into ``self.weyl_group``.
    """

```

```python

# DON'T: it is not a weight, or weight space, but a Weyl group
self._W_weight = self._finite_weight_space.weyl_group()
# DO
self.finite_weyl_group = self._finite_weight_space.weyl_group()
# DO: a bit verbose, use only to distinguish multiple Weyl groups in the same context
self.finite_weyl_group_from_weight_space = self._finite_weight_space.weyl_group()
```

```python
# DON'T: what is "ambient" space/set? "candidates" for what?
ambient_candidates = self.kl.affine_weyl_elements_from_bounded_translations(
            self.algebra,
            translations=normalized_translations,
            factor_order="st",
        )

# DO: if the logic involves only one set of `affine_weyl_elements`
affine_weyl_elements = self.kl.affine_weyl_elements_from_bounded_translations(
            self.algebra,
            translations=normalized_translations,
            factor_order="st",
        )
# DO: if the logic involves multiple sets of `affine_weyl_elements`, then add descriptive suffixes to distinguish them
affine_weyl_elements_from_bounded_translations = self.kl.affine_weyl_elements_from_bounded_translations(
            self.algebra,
            translations=normalized_translations,
            factor_order="st",
        )

# DO: use standard mathematical notation for sets
W_hat = self.kl.affine_weyl_elements_from_bounded_translations(
            self.algebra,
            translations=normalized_translations,
            factor_order="st",
        )
```