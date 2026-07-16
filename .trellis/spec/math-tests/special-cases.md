# Special Cases

List the small-rank, low-degree, degenerate, or exactly solvable cases that every implementation should pass.

For each case, include:

- input data
- expected output or property
- why this case matters mathematically

## Affine Verma Partition Cases

`pyw/tests/test_affine_verma_partition.py` must cover:

- $A_1^{(1)}$ low-coordinate coefficients, including $P(0)=1$;
- the $D_4^{(1)}$ coefficient of $\delta$, where finite-root partitions and
  four imaginary-root colors both contribute;
- invariance of lower coefficients when the coordinate box is enlarged;
- all 238 KL-support differences for the $D_4^{(1)}$ level $-2$, order-4
  benchmark, compared with an independently constructed recursive partition.

## Affine Freudenthal Cases

`pyw/tests/test_hybrid_freudenthal.py` must cover:

- an $A_1^{(1)}$ non-vacuum finite-dominant highest weight at nonzero grade;
- rejection of nonintegral finite Dynkin labels;
- $A_1^{(1)}$, $A_2^{(1)}$, and $B_2^{(1)}$ depth-one root domains;
- $B_2$ root/weight round trips, symmetric scalar products, and long/short
  root norms;
- exact shifted-norm detection of known degenerate and nondegenerate targets;
- explicit callback routing for every zero denominator;
- rejection of fractional callback results;
- the $D_4^{(1)}$ level $-2$ order-4 weight counts
  `1, 25, 169, 625, 1681`;
- exactly seven non-highest degenerate dominant targets;
- grade dimensions `1, 28, 329, 2632, 16380` after Freudenthal recursion.

## Hybrid Freudenthal-KL Cases

`pyw/tests/test_hybrid_freudenthal_kl.py` must cover the complete
$D_4^{(1)}$ level $-2$, order-4 calculation from an isolated KL cache:

- seven non-highest Freudenthal-degenerate dominant targets;
- 238 root-cone and nonzero-Verma requests;
- 81 distinct Q-tilde second arguments and 157 same-calculation cache hits;
- complete grade dimensions `1, 28, 329, 2632, 16380`.
