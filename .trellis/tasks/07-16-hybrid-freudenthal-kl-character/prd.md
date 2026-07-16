# Hybrid Freudenthal-KL Affine Character

## Goal

Add a fast affine-character computation path to `pyw` for supported affine highest-weight modules. The new path must reuse the existing order-bounded affine Weyl enumeration as a complete bounded database of KL orbit weights, but it must not compute every quotient inverse Kazhdan-Lusztig coefficient eagerly.

For a requested character order, the algorithm will:

1. enumerate exactly the affine Weyl representatives already selected by the current KL character implementation;
2. cache the correspondence between each representative and its shifted orbit weight;
3. compute ordinary weight multiplicities with Freudenthal recursion whenever its denominator is nonzero;
4. use KL coefficient extraction only at Freudenthal-degenerate target weights;
5. lazily compute and cache only the Qtilde values whose orbit weights can contribute to those targets;
6. assemble the complete character through the requested order.

The primary benchmark, used only for validation and performance measurement, is

$$
D_4^{(1)},\qquad \hat\lambda=-2\hat\omega_0,\qquad \operatorname{order}=4.
$$

## Background

The current `KazhdanLusztigCharacter.character_weight_list()` path constructs the complete bounded KL numerator:

$$
N_\lambda
=
\sum_w
\widetilde Q_{[w_\lambda],[w]}(1)e^{w\cdot\Lambda}.
$$

For the D4 vacuum benchmark at order 4, the current pipeline produces:

| Stage | Count |
|---|---:|
| Translations | 23 |
| Raw finite-Weyl-by-translation candidates | 4416 |
| Quotient representatives | 2208 |
| Representatives after Bruhat filtering | 2112 |
| Eager Qtilde evaluations requested | 2112 |

The affine Weyl enumeration itself is not the performance problem. Measured timings are approximately:

| Stage | Time |
|---|---:|
| Translation enumeration | 0.522 s |
| Translation-to-coroot conversion | 0.008 s |
| Finite Weyl by translation construction | 0.192 s |
| Strict affine Weyl enumeration total | 0.722 s |
| All preparation before Qtilde | 28.877 s |

The expensive stage is the eager computation of 2112 Qtilde values. A full `pyw` order-4 run did not finish within 20 minutes.

`scripts/negative_level.py` instead computes module weight multiplicities directly. For a target weight $\mu$, it uses Freudenthal recursion when

$$
F(\mu)
=
\lVert\lambda+\rho\rVert^2
-
\lVert\mu+\rho\rVert^2
\neq 0,
$$

and uses KL coefficient extraction only when $F(\mu)=0$.

For this benchmark only, after finite Weyl orbit reduction, the reference path has seven Freudenthal-degenerate target weights. Their KL-support sizes are

$$
12,\ 12,\ 12,\ 28,\ 58,\ 58,\ 58.
$$

This gives 238 Qtilde requests before deduplication and 81 distinct second arguments after deduplication.

## Mathematical Definitions

### KL Orbit Weight

For a bounded affine Weyl representative $w$, define

$$
\nu_w=w(\Lambda+\rho)-\rho.
$$

Every such weight satisfies

$$
\lVert\nu_w+\rho\rVert^2
=
\lVert\Lambda+\rho\rVert^2.
$$

These weights form the bounded KL numerator support. They are not the complete set of weights of the irreducible module.

### Freudenthal-Degenerate Target

A module weight $\mu$ is Freudenthal-degenerate when

$$
F(\mu)=0.
$$

This term must be used in code and documentation instead of ambiguously calling every such target a primitive null vector.

### KL Coefficient Extraction

At a Freudenthal-degenerate target $\mu$, its irreducible-module multiplicity is

$$
m(\mu)
=
\sum_w
\widetilde Q_{[w_\lambda],[w]}(1)
P(\nu_w-\mu),
$$

where $P$ is the affine Verma/Kostant partition multiplicity.

A representative can contribute only if

$$
\nu_w-\mu\in\widehat Q_+.
$$

This root-cone condition is the primary lazy-Qtilde query filter.

### Freudenthal Recursion

For $F(\mu)\neq0$,

$$
m(\mu)
=
\frac{
2\displaystyle\sum_{\alpha\in\widehat\Delta_+}
\operatorname{mult}(\alpha)
\sum_{j\geq1}
(\mu+j\alpha,\alpha)m(\mu+j\alpha)
}{F(\mu)}.
$$

All arithmetic must be exact. A computed multiplicity must be a nonnegative integer.

## Scope

### In Scope

- Reuse the existing order-bounded affine Weyl enumeration.
- Extract an orbit-preparation API that stops before Qtilde evaluation.
- Build exact bidirectional indexes between bounded representatives and orbit weights.
- Build the bounded module weight tree needed for Freudenthal recursion.
- Implement exact affine simple-root coordinates and root-cone membership checks.
- Implement bounded positive real and imaginary root enumeration.
- Implement affine Freudenthal recursion.
- Implement bounded affine Verma/Kostant partition multiplicities.
- Implement lazy KL coefficient extraction for Freudenthal-degenerate targets.
- Cache Qtilde by canonical Weyl-element key across all targets in one character computation.
- Assemble affine-character coefficients and, when supported by the finite Lie algebra backend, finite-part irreducible decompositions through the requested order.
- Preserve the existing eager KL numerator and character APIs unchanged.
- Validate the generic algorithm against the D4 level -2 vacuum benchmark through order 4.

### Out of Scope

- Replacing or removing the existing full KL numerator implementation.
- Claiming that uncomputed Qtilde values vanish.
- Recovering every omitted Qtilde value from the final character.
- Support for modules that are not highest-weight modules.
- Finite-Weyl canonicalization outside the supported finite-dominant-integral highest-weight domain.
- Support for twisted affine types in the first implementation.
- Optimizing orbit-index queries beyond a linear scan unless profiling shows it is necessary.
- Changing the mathematical convention for Qtilde or its evaluation at $q=1$.

## Technical Design

### 1. Bounded KL Orbit Preparation

Refactor the non-Qtilde portion of `KazhdanLusztigCharacter.character_weight_list()` into a reusable preparation method.

Proposed result object:

```python
@dataclass(frozen=True)
class BoundedKLOrbit:
    lambda_hat: AffineWeight
    Lambda_hat: AffineWeight
    w_to_lambda: Any
    stabilizer: tuple[Any, ...]
    representatives: tuple[Any, ...]
    weyl_to_weight: Mapping[WeylKey, AffineWeight]
    weight_to_weyl: Mapping[WeightKey, Any]
    by_grade: Mapping[int, tuple[WeightKey, ...]]
```

The preparation path must preserve the current semantics exactly:

1. find $\Lambda$ and $w_\lambda$;
2. use `translation_order = order_base + order`;
3. construct $W_{\mathrm{fin}}\times T$;
4. compute the full stabilizer;
5. compute shifted orbit weights;
6. deduplicate equal orbit weights, retaining the shortest representative;
7. apply the existing Bruhat condition;
8. return without calling Qtilde.

### 2. Exact Canonical Keys

Do not use floating-point data in new cache keys.

Proposed weight key:

```python
WeightKey = tuple[tuple[QQ, ...], QQ, QQ]
```

It contains exact affine Dynkin labels, level, and grade.

Proposed Weyl key:

```python
WeylKey = tuple[int, ...]
```

It contains a canonical reduced word using the same convention as the existing KL caches.

### 3. Affine Root-Cone Coordinates

Add a helper that solves

$$
\beta=\sum_{i=0}^r n_i\alpha_i
$$

for an affine weight difference $\beta$.

Required interface:

```python
def affine_simple_root_coordinates(
    algebra: AffineLieAlgebra,
    difference: AffineWeight,
) -> tuple[QQ, ...]:
    ...
```

Root-cone membership requires every coordinate to be an integer greater than or equal to zero.

### 4. Lazy KL Candidate Query

For a Freudenthal-degenerate target $\mu$, scan the bounded orbit representatives and retain only entries satisfying

$$
\nu_w-\mu\in\widehat Q_+.
$$

The first implementation may scan all bounded orbit weights for each target. For the D4 order-4 benchmark this costs at most

$$
7\times2112=14784
$$

exact root-cone tests, which is negligible compared with Qtilde evaluation.

The query order must be:

1. grade feasibility check;
2. affine root-cone check;
3. Verma partition lookup;
4. skip if the partition multiplicity is zero;
5. lazily compute Qtilde only after all cheap filters pass.

### 5. Lazy Qtilde Cache

Use one cache for the complete character calculation:

```python
q_tilde_cache: dict[WeylKey, Any]
```

The fixed first argument and stabilizer belong to the surrounding character context. If the cache is persisted or shared between contexts, its key must additionally include the first Weyl argument, Cartan type, and stabilizer identity.

### 6. Bounded Weight Tree

Generate module weights through the requested order and group them by relative affine depth. For a highest weight that is dominant integral along the finite simple roots, use finite-simple-root lowering closure at each depth and use $\alpha_0$ lowering to advance the affine depth.

Use sets for deduplication. Preserve exact affine weights.

The calculation order must follow the positive-root partial order, implemented through the height of

$$
\lambda-\mu=\sum_i n_i\alpha_i.
$$

Higher weights must be computed before lower weights.

### 7. Bounded Positive Roots

Generate the positive roots needed through the requested grade:

$$
\alpha+n\delta,
\qquad
-\alpha+n\delta,
\qquad
n\delta.
$$

Real-root multiplicity is 1. The multiplicity of $n\delta$ is the rank of the finite Lie algebra.

Do not depend on the Sage 10.4 extended-weight-lattice `positive_real_roots()` iterator that currently raises a parent-coercion error.

### 8. Affine Verma Partition Table

Build one bounded partition table per character computation and reuse it for all KL targets:

```python
class BoundedAffineVermaPartition:
    def multiplicity(self, coordinates: tuple[int, ...]) -> int:
        ...
```

The table must include real and imaginary root multiplicities and must be truncated by the maximum grade and root height actually required by the target set.

### 9. Hybrid Multiplicity Engine

Proposed class:

```python
class KazhdanLusztigFreudenthalCharacter:
    def multiplicities(
        self,
        lambda_hat: AffineWeight,
        *,
        order: int,
    ) -> dict[int, dict[AffineWeight, int]]:
        ...

    def character(
        self,
        lambda_hat: AffineWeight,
        *,
        order: int,
    ) -> Any:
        ...
```

Keep this implementation in a dedicated module or a clearly separated class. Do not add the full algorithm directly into `kazhdan_lusztig.py`, whose responsibility is polynomial computation.

### 10. Main Algorithm

```python
orbit = prepare_bounded_kl_orbit(lambda_hat, order=order)
weights = build_bounded_weight_tree(lambda_hat, order=order)
partition = build_bounded_verma_partition(weights, orbit)

multiplicities = {weight_key(lambda_hat): 1}

for mu in targets_sorted_by_root_height(weights):
    if mu == lambda_hat:
        continue

    denominator = shifted_norm(lambda_hat) - shifted_norm(mu)

    if denominator != 0:
        value = freudenthal_multiplicity(mu, multiplicities, weights)
    else:
        value = kl_coefficient_at_weight(mu, orbit, partition)

    require_nonnegative_integer(value)
    multiplicities[weight_key(mu)] = value

return assemble_character(multiplicities)
```

## Implementation Plan

### Phase 1: Extract the Bounded Orbit API

- Extract current character steps 1 through 7 into a reusable preparation method.
- Add `BoundedKLOrbit` and exact canonical keys.
- Preserve `character_weight_list()` behavior by making it consume the prepared orbit and then perform eager Qtilde evaluation as before.
- Add instrumentation for counts and timing.

Deliverable: a no-Qtilde bounded orbit with exactly the same representatives as the current eager path.

### Phase 2: Root Coordinates and Candidate Queries

- Implement affine simple-root coordinate conversion.
- Implement exact $\widehat Q_+$ membership.
- Implement `candidates_above(mu)`.
- Compare candidate lists against `scripts/negative_level.py.find_nulls_above()`.

Benchmark deliverable: the seven D4 test targets produce support sizes $12,12,12,28,58,58,58$, with 81 unique representatives.

### Phase 3: Bounded Roots and Verma Multiplicities

- Implement bounded positive real and imaginary roots.
- Implement the bounded affine Verma partition table.
- Validate partition multiplicities against the reference script on all KL-support differences in the D4 order-4 benchmark.

Deliverable: exact Verma multiplicities for every queried $\nu_w-\mu$.

### Phase 4: Freudenthal Multiplicity Engine

- Implement the bounded weight domain for supported affine highest-weight modules whose finite Dynkin labels are dominant integral.
- Implement target ordering by affine root height.
- Implement Freudenthal recursion with exact arithmetic.
- Detect denominator-zero targets and route them to the KL branch.
- Reject nonintegral or negative outputs.

Benchmark deliverable: all nondegenerate D4 test-target multiplicities match the reference through order 4.

### Phase 5: Hybrid Character Assembly

- Assemble finite Weyl orbit multiplicities by grade.
- Decompose each grade into irreducible characters of the corresponding finite Lie algebra whenever that decomposition backend is available.
- Expose the hybrid public API without changing the existing eager API.
- Add profiling counters for target counts, query hits, unique Qtilde evaluations, cache hits, and phase timing.

Benchmark deliverable: the generic API reproduces the complete D4 level -2 vacuum character through order 4.

### Phase 6: Performance and Documentation

- Benchmark fresh and cached runs.
- Document mathematical terminology and algorithm limitations.
- Update implementation, definition, and test specs with executable contracts.

## Acceptance Criteria

### Functional

- [ ] A bounded-orbit API returns the representatives selected by the current strict `order`-bounded KL path without computing Qtilde.
- [ ] Existing eager `KazhdanLusztigCharacter` behavior and output remain unchanged.
- [ ] Exact bidirectional representative/weight indexes are available.
- [ ] Freudenthal recursion is used only when its denominator is nonzero.
- [ ] Every denominator-zero target is routed to KL coefficient extraction.
- [ ] KL coefficient extraction filters by $\nu_w-\mu\in\widehat Q_+$ before Qtilde evaluation.
- [ ] Verma partition multiplicity is checked before Qtilde evaluation.
- [ ] Qtilde values are computed lazily and cached across target weights.
- [ ] All returned multiplicities are exact nonnegative integers.
- [ ] The new API returns the complete character through the requested order, not the full explicit KL numerator.

### D4 Mathematical Benchmark — Test Data Only

- [ ] For $D_4^{(1)}$, $\hat\lambda=-2\hat\omega_0$, order 4, strict enumeration gives 23 translations, 4416 raw candidates, 2208 quotient representatives, and 2112 Bruhat-valid representatives.
- [ ] Orbit preparation performs zero Qtilde evaluations.
- [ ] Exactly seven finite-Weyl target representatives are Freudenthal-degenerate through order 4.
- [ ] Their candidate support sizes are $12,12,12,28,58,58,58$.
- [ ] The union of the candidate Weyl representatives contains 81 distinct second arguments before any stronger partition-zero filtering.
- [ ] Grade dimensions are

  $$
  1,\ 28,\ 329,\ 2632,\ 16380.
  $$

- [ ] The graded dimension series is

  $$
  1+28q+329q^2+2632q^3+16380q^4+O(q^5).
  $$

- [ ] The expected finite D4 decomposition for this benchmark is:

  $$
  \begin{aligned}
  q^0:&\quad \mathbf 1,\\
  q^1:&\quad V(\omega_2),\\
  q^2:&\quad V(2\omega_2)+V(\omega_2)+\mathbf 1,\\
  q^3:&\quad V(3\omega_2)+V(2\omega_2)+V(\omega_1+\omega_3+\omega_4)+2V(\omega_2)+\mathbf 1,\\
  q^4:&\quad V(4\omega_2)+V(3\omega_2)+V(\omega_1+\omega_2+\omega_3+\omega_4)\\
  &\qquad +3V(2\omega_2)+V(\omega_1+\omega_3+\omega_4)\\
  &\qquad +V(2\omega_1)+V(2\omega_3)+V(2\omega_4)+3V(\omega_2)+2\mathbf 1.
  \end{aligned}
  $$

### Performance

- [ ] For the D4 order-4 benchmark, the hybrid path evaluates no more than 81 distinct Qtilde second arguments before optional stronger filtering.
- [ ] The implementation reports Qtilde request count, unique evaluation count, cache-hit count, and phase timings.
- [ ] Fresh D4 order-4 runtime target: less than 60 seconds on the benchmark machine.
- [ ] Cached D4 order-4 runtime target: less than 10 seconds on the benchmark machine.
- [ ] Strict affine Weyl enumeration remains below 2 seconds on the benchmark machine.
- [ ] No test or benchmark relies on preexisting user KL cache data unless explicitly testing warm-cache behavior.

### Quality

- [ ] New mathematical terms are documented and used consistently.
- [ ] Cache keys use exact values and canonical Weyl words.
- [ ] The implementation does not silently map a general nonintegrable weight to a finite dominant representative.
- [ ] Inputs outside the documented applicability domain fail with a clear error instead of producing a plausible but unverified affine character.
- [ ] Existing KL, affine-weight, and affine-Weyl tests pass.

## Test Plan

### Unit Tests: Exact Weight and Root Operations

Add tests for:

- affine weight canonical keys;
- Weyl reduced-word keys;
- affine simple-root coordinate reconstruction;
- root-cone membership with positive, boundary-zero, fractional, and negative coordinates;
- shifted norm and Freudenthal denominator using exact `QQ` arithmetic;
- bounded positive real and imaginary root generation;
- real and imaginary root multiplicities.

### Unit Tests: Bounded Orbit Preparation

Test that:

- preparation calls no Qtilde method;
- prepared representatives match the current eager path exactly by reduced word;
- prepared orbit weights match exactly by weight key;
- duplicate orbit weights retain the shortest representative;
- the identity is present in the stabilizer;
- the `order_base` shift is applied correctly;
- the D4 order-4 benchmark produces the expected counts.

### Unit Tests: Lazy Candidate Query

For fixed target weights, test:

- every returned $\nu_w-\mu$ lies in $\widehat Q_+$;
- every omitted bounded orbit weight fails the root-cone condition;
- boundary cases with one or more zero simple-root coefficients are retained;
- query results are independent of iteration order;
- repeated targets reuse the same canonical Qtilde cache entries.

### Unit Tests: Verma Partition

Test:

- $P(0)=1$;
- differences outside $\widehat Q_+$ have multiplicity zero;
- low-height real-root differences have the expected multiplicity;
- imaginary-root multiplicities use the finite rank;
- all 238 queries in the D4 benchmark support set match the reference `negative_level.py` Verma multiplicities;
- rebuilding with a larger bound preserves all lower-bound coefficients.

### Unit Tests: Freudenthal

Test:

- highest-weight multiplicity is 1;
- weights outside the highest-weight cone return zero or fail validation as specified;
- nondegenerate targets match known reference multiplicities;
- denominator-zero targets never divide and always invoke the KL branch;
- every result is an exact nonnegative integer;
- dependency ordering computes every required higher weight first;
- recursion-loop detection reports an implementation error.

### Integration Tests: Hybrid KL Branch

For the seven degenerate targets in the D4 benchmark:

- compare candidate support sets with the reference script;
- compare individual Qtilde values for all 81 unique representatives;
- compare Verma-weighted KL sums;
- compare final target multiplicities;
- verify that duplicate support representatives cause cache hits rather than recomputation.

### Integration Tests: Complete Character

Run the generic affine-character API on the D4 level -2 vacuum benchmark at orders 0, 1, 2, 3, and 4. For each order:

- compare all computed dominant weight multiplicities with the reference script;
- compare finite irreducible decompositions;
- compare grade dimensions;
- verify truncation consistency: order $n+1$ must preserve every coefficient through order $n$.

### Existing-Path Compatibility Tests

- Run the existing affine KL character tests.
- Compare eager-path prepared representatives before and after refactoring.
- Confirm that direct eager Qtilde output is unchanged on small A2 cases.
- Confirm persistent Q and Qtilde cache formats remain compatible unless an explicit migration is added.

### Performance Tests

Use isolated temporary cache directories.

Measure:

1. cold bounded-orbit preparation;
2. cold hybrid character;
3. same-process warm hybrid character;
4. new-process persistent-cache hybrid character;
5. number of candidate-query hits;
6. number of unique Qtilde evaluations;
7. Qtilde cache hits;
8. Freudenthal target count;
9. phase-level wall time.

Performance tests should report thresholds but remain opt-in if the shared test environment is too variable for strict wall-clock assertions.

### Required Commands

Core tests must use Sage:

```bash
sage -python -m pytest pyw/tests/test_affine_kl_context.py pyw/tests/test_affine_kl_character.py
```

Add focused commands for the new test files, for example:

```bash
sage -python -m pytest pyw/tests/test_hybrid_affine_character.py
sage -python -m pytest pyw/tests/test_affine_verma_partition.py
```

## Risks and Safeguards

### Mathematical Terminology

Risk: conflating Freudenthal-degenerate weights, primitive singular-vector weights, and KL orbit weights.

Safeguard: use distinct names and document each definition.

### Incomplete Bounded Orbit

Risk: a denominator-zero target has a root-cone-supported orbit weight that was omitted by an incorrect order bound.

Safeguard: preserve the existing strict enumeration semantics and add a coverage assertion against reference support sets.

### Incorrect Finite Weyl Reduction

Risk: applying finite Weyl orbit reduction where the module character is not finite-Weyl invariant.

Safeguard: require an untwisted affine highest-weight input whose finite Dynkin labels are dominant integral; reject inputs outside this explicitly documented applicability domain.

### Sage Root Iterator Compatibility

Risk: Sage 10.4 cannot iterate the extended-weight-lattice positive real roots used by the reference script.

Safeguard: construct bounded affine roots explicitly from finite roots and $\delta$.

### Qtilde Convention Drift

Risk: polynomial-first evaluation and direct evaluation at $q=1$ disagree because of backend behavior.

Safeguard: preserve the current project convention, add representative-level comparisons, and do not mix cache namespaces.

### Affine Character Versus Numerator API

Risk: callers assume the hybrid path returns the complete explicit KL numerator rather than an affine character.

Safeguard: use a distinct class/API and document that the output is the affine character through a requested order.

## Related Files

- `pyw/core/character.py`
- `pyw/core/kazhdan_lusztig.py`
- `pyw/core/affine_weight.py`
- `pyw/core/affine_lie_algebra.py`
- `pyw/core/hybrid_affine_character.py`
- `pyw/core/weyl_group.py`
- `scripts/negative_level.py`
- `pyw/tests/test_affine_kl_context.py`
- `pyw/tests/test_affine_kl_character.py`
- `pyw/tests/test_hybrid_affine_character.py`

## Definition of Done

- The implementation and all required tests satisfy the acceptance criteria.
- The generic implementation reproduces the D4 order-4 benchmark from a clean cache environment.
- Cold and warm benchmark reports are recorded in the task directory.
- Existing eager KL behavior remains available and tested.
- Relevant definition, implementation, core-test, and math-test specs are updated with executable contracts.
- The exact applicability domain and unsupported inputs are documented.

## Implementation Progress

### 2026-07-16: Phases 1 and 2 Complete

Implemented the reusable bounded KL orbit preparation path in
`pyw/core/character.py`:

- `BoundedKLOrbit` stores the original and dominant weights, fixed first Weyl
  argument, stabilizer, strict Bruhat-valid representatives, exact
  representative-to-weight and weight-to-representative indexes, and a grade
  index.
- `KazhdanLusztigCharacter.prepare_bounded_kl_orbit()` executes the existing
  eager path through Bruhat filtering without evaluating Qtilde.
- `character_weight_list()` consumes the prepared orbit and preserves the
  existing Qtilde call order, stabilizer argument, zero-coefficient filtering,
  and output format.

Implemented exact root-cone and candidate-query support in
`pyw/core/hybrid_affine_character.py`:

- `affine_simple_root_coordinates()` uses
  $n_0=\operatorname{grade}(\beta)$ and $n_i=c_i+a_i n_0$ for
  $\alpha_0=\delta-\theta$.
- `is_in_affine_positive_root_cone()` requires exact nonnegative integral
  coordinates and level zero.
- `BoundedKLOrbit.candidates_above()` applies the grade feasibility check and
  affine root-cone filter before later Verma and Qtilde work.
- The initial implementation explicitly rejects twisted affine types.

Added exact `AffineWeightKey` support and fixed two conversion issues discovered
during Codex review and evidence-based debugging:

- affine Dynkin label reconstruction now uses comarks for the level relation;
- `AffineLieAlgebra.from_sagemath()` reads an appended grade only from extended
  affine weight vectors, not ordinary affine Dynkin-label vectors.

### D4 Order-4 Benchmark Verification

For $D_4^{(1)}$, $\hat\lambda=-2\hat\omega_0$, and order $4$:

| Measurement | Result |
|---|---:|
| Bruhat-valid representatives | 2112 |
| Seven candidate support sizes | 12, 12, 12, 28, 58, 58, 58 |
| Distinct candidate Weyl keys | 81 |
| Qtilde evaluations during orbit preparation | 0 |

Within this benchmark, the seven non-highest Freudenthal-degenerate dominant targets are the triality
triples at grades 2 and 4 together with
$\omega_1+\omega_3+\omega_4$ at grade 3.

### Verification

```bash
sage -python -m pytest --no-cov \
  pyw/tests/test_hybrid_affine_character.py \
  pyw/tests/test_affine_kl_character.py \
  pyw/tests/test_affine_weight.py \
  pyw/tests/test_affine_lie_algebra.py
```

Result: **108 passed**.

The final focused phase-two verification produced **15 passed**, including the
D4 order-4 support benchmark. Codex reviewed each completed step; all reported
correctness issues were fixed and re-tested.

### 2026-07-16: Phase 3 Complete

Implemented bounded untwisted affine positive roots and affine Verma/Kostant
partition multiplicities in `pyw/core/hybrid_affine_character.py`:

- `bounded_affine_positive_roots(algebra, max_grade)` explicitly constructs
  the real roots $\alpha+n\delta$ and $-\alpha+n\delta$, and the imaginary
  roots $n\delta$, without using Sage's extended affine root iterator;
- `BoundedAffinePositiveRoot` stores exact affine simple-root coordinates,
  root multiplicity, and whether the root is imaginary;
- real roots have multiplicity one and each $n\delta$ has multiplicity equal
  to the finite rank;
- `BoundedAffineVermaPartition(algebra, max_coordinates)` builds a sparse exact
  partition table inside a componentwise affine-simple-root coordinate box;
- nonintegral bounds, negative bounds, finite Cartan types, and twisted affine
  types fail explicitly instead of producing truncated or unverified data.

For the 238 support differences in the D4 benchmark, the required coordinate box is
`(4, 2, 4, 2, 2)`. It contains 51 bounded positive-root records and produces a
sparse table with 675 nonzero entries. Every one of the 238 queried
multiplicities agrees with an independently constructed recursive partition
calculation.

Verification:

```bash
sage -python -m pytest --no-cov pyw/tests/test_affine_verma_partition.py
sage -python -m pytest --no-cov \
  pyw/tests/test_hybrid_affine_character.py \
  pyw/tests/test_affine_kl_character.py
ruff check \
  pyw/core/hybrid_affine_character.py \
  pyw/core/__init__.py \
  pyw/tests/test_affine_verma_partition.py
```

Results: **6 passed**, **15 passed**, and **all Ruff checks passed**.

### 2026-07-16: Phase 4 Complete

Implemented the bounded affine highest-weight domain and exact Freudenthal engine in
`pyw/core/hybrid_affine_character.py`:

- `build_bounded_weight_tree()` advances depth by subtracting $\alpha_0$ and
  closes each grade under finite-simple-root lowering;
- exact `AffineWeightKey` values provide all deduplication and lookup keys;
- targets are processed in increasing affine-root height of
  $\lambda-\mu$, so every higher-weight dependency is already available;
- `classically_dominant_affine_weight()` maps each finite Weyl orbit to its
  dominant representative while preserving level and grade;
- `freudenthal_denominator()` uses the exact shifted affine norm with
  $\widehat\rho$;
- `compute_bounded_freudenthal_multiplicities()` applies real and imaginary
  root multiplicities, routes every zero denominator through an explicit
  callback, and rejects fractional or negative results.

The bounded weight-domain construction accepts arbitrary affine level and
highest-weight grade. Its module-theoretic condition is finite-simple-root
dominant integrality, which is exactly the condition needed for finite Weyl
orbit reduction. Fractional or negative finite Dynkin labels fail explicitly.

The generality review also corrected two non-simply-laced foundations:

- root-to-weight conversion uses Cartan-matrix columns, with the inverse
  conversion using the matching inverse-matrix rows;
- the finite invariant form uses the root-length symmetrizer rather than the
  generally nonsymmetric matrix $C^{-1}$ alone.

These corrections are verified on $B_2$ simple roots, highest roots,
fundamental weights, and the depth-one affine weight domain. The D4 data below
remain benchmark assertions only and do not enter production control flow.

As a benchmark-only verification, for the D4 level $-2$ vacuum through order 4:

| Measurement | Result |
|---|---:|
| Weight counts by grade | 1, 25, 169, 625, 1681 |
| Dominant targets by grade | 1, 2, 7, 15, 30 |
| Non-highest zero-denominator targets | 7 |
| Grade dimensions | 1, 28, 329, 2632, 16380 |

For this benchmark test, the seven callback values were independently extracted from the finite D4
decomposition in this PRD: `2, 2, 2, 8, 16, 16, 16`. All remaining dominant
weight multiplicities were produced by Freudenthal recursion.

Verification:

```bash
sage -python -m pytest --no-cov \
  pyw/tests/test_hybrid_freudenthal.py \
  pyw/tests/test_hybrid_affine_character.py \
  pyw/tests/test_affine_kl_character.py
mypy pyw/core/hybrid_affine_character.py pyw/core/__init__.py
ruff check \
  pyw/core/hybrid_affine_character.py \
  pyw/core/__init__.py \
  pyw/tests/test_hybrid_freudenthal.py
```

Results: **21 passed**, **no mypy issues**, and **all Ruff checks passed**.

### 2026-07-16: Generality Review

Reviewed Phases 1 through 4 to remove assumptions tied to a particular Lie
algebra or highest weight:

- removed the $k\widehat\Lambda_0$ and grade-zero input restriction;
- required finite-simple-root dominant integrality instead;
- added exact B2 root/weight and invariant-form support;
- distinguished bounded-domain membership from an uncomputed Freudenthal
  dependency;
- made explicit `Lambda_hat` and `w_to_lambda` inputs an all-or-nothing pair;
- changed dominant-Lambda checks to use all affine Dynkin labels instead of a
  positional vector prefix;
- retained D4 level $-2$ and seven-target data exclusively as benchmark fixtures and assertions.

Verification:

```bash
sage -python -m pytest --no-cov \
  pyw/tests/test_hybrid_affine_generality.py \
  pyw/tests/test_hybrid_freudenthal.py \
  pyw/tests/test_hybrid_affine_character.py \
  pyw/tests/test_affine_verma_partition.py \
  pyw/tests/test_affine_kl_character.py \
  pyw/tests/test_affine_lie_algebra.py \
  pyw/tests/test_affine_weight.py
```

Result: **126 passed**. A final focused run after the bounded-orbit API changes
produced **34 passed**. Mypy and Ruff checks for the new hybrid modules pass.

### Current Work

Phase 5 is next:

- connect the zero-denominator callback to bounded KL candidates;
- filter candidates by affine Verma multiplicity;
- lazily compute and cache Qtilde values;
- expose `KazhdanLusztigFreudenthalCharacter` and assemble the complete character.

### 2026-07-16: Phase 5 Complete

Implemented `KazhdanLusztigFreudenthalCharacter` and the complete lazy KL branch in
`pyw/core/hybrid_affine_character.py`:

- degenerate targets are filtered by grade, affine root cone, and bounded
  Verma multiplicity before Qtilde evaluation;
- one exact reduced-word cache is shared across all targets;
- prepared orbit entries cache representative keys and affine root coordinates;
- candidate support sets are computed once and reused by KL extraction;
- `multiplicities()` returns dominant multiplicities by relative depth;
- `character()` expands them to every bounded module weight;
- `profile_stats()` reports target, request, cache, and phase timing counters.

The isolated-cache $D_4^{(1)}$ level $-2$, order-4 run produced:

| Measurement | Result |
|---|---:|
| Grade dimensions | `1, 28, 329, 2632, 16380` |
| Degenerate targets | 7 |
| Candidate and Verma-nonzero requests | 238 |
| Unique Qtilde evaluations | 81 |
| Same-calculation Qtilde cache hits | 157 |
| Orbit preparation | 46.91 s |
| Partition preparation | 0.06 s |
| Hybrid multiplicity calculation | 1.33 s |
| Total cold runtime | 48.64 s |

This satisfies the cold runtime target of less than 60 seconds on the benchmark
machine. The focused orbit, eager compatibility, and hybrid tests produced
**19 passed** before the final order-4 assertion was promoted into the test
suite. Ruff and mypy pass for the new hybrid module and exports; the existing
`character.py` file retains unrelated historical lint and type findings.

### 2026-07-16: Phase 6 Benchmark Record

`HybridCharacterStats` now includes underlying Q-cache counters. The benchmark
script `benchmark_phase6.py` uses an isolated cache directory and confirms the
same D4 dimensions in every mode:

| Mode | Total time | Q calls | Q cache hits |
|---|---:|---:|---:|
| In-process cold | 47.80 s | 162 | 0 |
| In-process second run | 46.78 s | 162 | 162 |
| Isolated persistent-cache cold process | 81.86 s | 162 | 0 |
| Isolated persistent-cache second process | 68.47 s | 162 | 162 |

The Q cache is correctly reused both within and across processes. The complete
warm runtime does **not** meet the 10-second target because strict bounded
orbit preparation is rebuilt for every public call and dominates runtime. No
incorrect performance claim is made; persistent bounded-orbit reuse is deferred
until a separate API and invalidation contract are designed.

### 2026-07-16: Persistent Bounded-Orbit Cache

Implemented the deferred cache as an opt-in `orbit_cache_dir` argument on
`KazhdanLusztigCharacter.prepare_bounded_kl_orbit()` and
`KazhdanLusztigFreudenthalCharacter`. It stores only reconstructible data in a
SQLite/WAL database: exact affine Dynkin labels and grade, canonical reduced
words, stabilizer words, and relative affine-root coordinates. It never stores
Sage objects or Qtilde values. Caller-provided translations are rejected because
they do not provide a stable cache identity.

The D4 order-4 snapshot round trip reconstructs all 2112 representatives in
**5.56 s**, replacing the approximately 47-second strict-orbit rebuild. The
new focused cache test verifies an A2 round trip, exact reconstructed words and
weights, and that a hit skips translation enumeration.

### 2026-07-16: Complete Warm Benchmark

With an isolated directory containing both persistent Q values and a bounded
orbit snapshot, a new D4 order-4 process produced the complete character in
**5.49 s**, meeting the less-than-10-second warm target. The cold process took
**50.47 s**. Both produced grade dimensions
`1, 28, 329, 2632, 16380`; the warm process reconstructed 2112 orbit
representatives in **4.14 s** and had **162/162** underlying Q calls served by
the persistent cache.
