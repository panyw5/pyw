"""Quantum Drinfeld-Sokolov reduction workflow (MVP for :math:`sl_n`).

Core scope in this module:
1) Build reduction data from a partition-labelled nilpotent orbit.
2) Enumerate a first-pass principal-admissible highest-weight pool.
3) Apply generic linear non-degeneracy constraints on affine Dynkin labels.
4) Compute central charge and reduced highest-weight conformal weight.
5) Evaluate characters for MVP-supported cases by reusing theta utilities.

Concrete case studies should live in ``demos/`` (not in this algorithm module).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from sage.all import QQ, CartanType, LieAlgebra, Matrix

from pyw.algorithms.nilpotent_orbits import NilpotentOrbit, nilpotent_orbits
from pyw.embedding.sl2_triple import RootSpaceGrading, SL2Triple
from pyw.fractional.level import FractionalLevel
from pyw.utils.theta_functions import sl2_boundary_character


def _to_fraction(value: Any) -> Fraction:
    """Convert a Sage/Python rational-like value to :class:`fractions.Fraction`."""
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
        num = value.numerator() if callable(value.numerator) else value.numerator
        den = value.denominator() if callable(value.denominator) else value.denominator
        return Fraction(int(num), int(den))
    return Fraction(value)


@dataclass(frozen=True)
class PrincipalAdmissibleLevelInput:
    """Input parameters for a principal-admissible level ``k=-h^vee+p/u``."""

    p: int
    u: int


@dataclass(frozen=True)
class ReductionData:
    """Nilpotent/triple/grading data used by DS reduction."""

    orbit: NilpotentOrbit
    triple: SL2Triple
    grading: RootSpaceGrading
    weighted_dynkin: Tuple[int, ...]
    h_coroot_coeffs: Tuple[Fraction, ...]
    x_coroot_coeffs: Tuple[Fraction, ...]
    h_f_dimension: int
    delta_l_size: int


@dataclass(frozen=True)
class LinearNonDegeneracyCondition:
    """Linear condition on affine Dynkin labels.

    For labels ``(lambda_0, ..., lambda_r)``, evaluates
    ``lhs = sum_i coeffs[i] * lambda_i`` and compares to ``rhs``.
    Supported relations: ``"neq"`` and ``"eq"``.
    """

    coeffs: Tuple[Fraction, ...]
    rhs: Fraction
    relation: str = "neq"


def is_nondegenerate_affine_dynkin(
    labels: Sequence[Fraction], conditions: Sequence[LinearNonDegeneracyCondition]
) -> bool:
    """Return ``True`` iff all non-degeneracy conditions are satisfied."""
    labels_f = tuple(_to_fraction(v) for v in labels)
    for cond in conditions:
        if len(cond.coeffs) != len(labels_f):
            raise ValueError("Condition dimension does not match affine label dimension")
        lhs = sum(cond.coeffs[i] * labels_f[i] for i in range(len(labels_f)))
        if cond.relation == "neq":
            if lhs == cond.rhs:
                return False
        elif cond.relation == "eq":
            if lhs != cond.rhs:
                return False
        else:
            raise ValueError(f"Unsupported relation: {cond.relation}")
    return True


@dataclass(frozen=True)
class AdmissibleWeightCandidate:
    """A finite-part highest-weight candidate in Dynkin-label form."""

    finite_dynkin_labels: Tuple[int, ...]
    lambda0_coeff: int
    affine_level: int
    reduction_grade: int
    cohomology_computed_upto_grade: int = 0
    cohomology_h0_dim: int = 0
    cohomology_h1_dim: int = 0
    cohomology_by_grade: Tuple[Tuple[int, int], ...] = ((0, 0),)
    cohomology_hq_by_grade: Tuple[Tuple[int, Tuple[int, int]], ...] = ((0, (0, 0)),)
    survives_reduction: bool = False


@dataclass(frozen=True)
class ReducedModuleData:
    """Computed reduced-module observables for a candidate weight."""

    candidate: AdmissibleWeightCandidate
    central_charge: Fraction
    conformal_weight: Fraction


@dataclass
class CharacterEvaluation:
    """Character evaluation result for one reduced module."""

    supported: bool
    value: Optional[Any] = None
    reason: str = ""


@dataclass
class QuantumDSReductionResult:
    """Structured output of the MVP DS reduction workflow."""

    cartan_type: str
    rank: int
    partition: Tuple[int, ...]
    level: FractionalLevel
    reduction_data: ReductionData
    all_candidates: List[AdmissibleWeightCandidate]
    surviving_candidates: List[AdmissibleWeightCandidate]
    all_modules: List[ReducedModuleData]
    surviving_modules: List[ReducedModuleData]

    def evaluate_character(
        self,
        module_index: int,
        tau: Any,
        z: Any,
        num_terms: int = 100,
    ) -> CharacterEvaluation:
        """Evaluate reduced character for MVP-supported cases.

        Supported path:
        - ``A_1`` + principal orbit + boundary-style level ``p=2`` + odd ``u``.
          Uses :func:`pyw.utils.theta_functions.sl2_boundary_character`.
        """
        if module_index < 0 or module_index >= len(self.all_modules):
            return CharacterEvaluation(False, reason="module_index out of range")

        module = self.all_modules[module_index]

        if self.cartan_type != "A" or self.rank != 1:
            return CharacterEvaluation(False, reason="MVP character only supports sl2 path")

        if tuple(self.partition) != (2,):
            return CharacterEvaluation(
                False,
                reason="MVP sl2 character path currently requires principal partition [2]",
            )

        if int(self.level.p) != 2:
            return CharacterEvaluation(False, reason="MVP sl2 character path requires boundary p=2")

        u = int(self.level.u)
        if u % 2 == 0:
            return CharacterEvaluation(False, reason="sl2 boundary character requires odd u")

        j = module.candidate.finite_dynkin_labels[0] % u
        value = sl2_boundary_character(u=u, j=j, tau=tau, z=z, num_terms=num_terms)
        return CharacterEvaluation(True, value=value)


def _dominant_finite_dynkin_labels(rank: int, affine_level: int) -> List[Tuple[int, ...]]:
    """Enumerate finite Dynkin labels ``(a_1,...,a_rank)`` with ``sum a_i <= level``."""
    labels: List[Tuple[int, ...]] = []
    for coeffs in product(range(affine_level + 1), repeat=rank):
        if sum(coeffs) <= affine_level:
            labels.append(tuple(int(c) for c in coeffs))
    return labels


def _quadratic_form_A_rank(rank: int) -> Matrix:
    """Inverse Cartan matrix for finite ``A_rank`` as a rational matrix."""
    cartan = Matrix(QQ, CartanType(["A", rank]).cartan_matrix())
    return cartan.inverse()


def _finite_weight_norm_A(labels: Sequence[int], rank: int) -> Fraction:
    """Compute ``(lambda, lambda+2rho)`` for finite ``A_rank`` from Dynkin labels."""
    a = Matrix(QQ, rank, 1, list(labels))
    ones = Matrix(QQ, rank, 1, [1] * rank)
    inv_c = _quadratic_form_A_rank(rank)
    value = (a.transpose() * inv_c * (a + 2 * ones))[0, 0]
    return _to_fraction(value)


def _build_reduction_data(cartan_type: str, rank: int, partition: Sequence[int]) -> ReductionData:
    """Build nilpotent/triple/grading data by reusing existing APIs."""
    partition_tuple = tuple(sorted((int(x) for x in partition), reverse=True))
    all_orbits = nilpotent_orbits(cartan_type, rank)
    orbit = next((o for o in all_orbits if tuple(o.partition) == partition_tuple), None)
    if orbit is None:
        raise ValueError(f"Partition {list(partition_tuple)} not valid for {cartan_type}_{rank}")

    triple = orbit.sl2_triple()
    grading = orbit.root_grading()
    h_coeffs = tuple(_to_fraction(c) for c in triple.h_coroot_coeffs)
    x_coeffs = tuple(c / 2 for c in h_coeffs)
    h_f_dim = len(triple.h_f())
    delta_l_size = len(triple.delta_l())

    return ReductionData(
        orbit=orbit,
        triple=triple,
        grading=grading,
        weighted_dynkin=tuple(int(d) for d in triple.wdd),
        h_coroot_coeffs=h_coeffs,
        x_coroot_coeffs=x_coeffs,
        h_f_dimension=h_f_dim,
        delta_l_size=delta_l_size,
    )


def _principal_central_charge_sl_n(n: int, k: Fraction) -> Fraction:
    """Principal DS central charge for ``W_k(sl_n)``."""
    n_frac = Fraction(n, 1)
    return (n_frac - 1) - Fraction(n * (n * n - 1), 1) * ((k + n_frac - 1) ** 2) / (k + n_frac)


def _compute_reduction_candidates(
    rank: int,
    affine_level: int,
    weighted_dynkin: Sequence[int],
    reduction_data: ReductionData,
    cohomology_max_grade: int = 0,
    nondegeneracy_conditions: Optional[Sequence[LinearNonDegeneracyCondition]] = None,
) -> List[AdmissibleWeightCandidate]:
    """Build candidates and compute phase-1 BRST cohomology at grade ``N=0``.

    Phase-1 policy computes a finite-dimensional BRST complex block
    ``C^0 -> C^1 -> C^2`` on the top-weight sector:

    - ``C^0``: span of highest-weight state (dimension 1)
    - ``C^1``: ghosts dual to roots in ``g_{>0}`` from the DS grading
    - ``C^2``: wedge products of two ``C^1`` ghosts

    Differential model:
    - ``d_0`` uses simple-root BRST constraints on highest weight,
      with coefficients ``lambda_i - chi_i`` where ``chi_i`` is the DS
      character support extracted from the nilpotent ``f``.
    - ``d_1`` is the Chevalley-Eilenberg differential on root-dual basis,
      using root-addition incidence (MVP normalization ``N_{beta,gamma}=1``).

    Cohomology dimensions are computed exactly by linear algebra:
    ``H^0 = ker(d_0)``, ``H^1 = ker(d_1)/im(d_0)``.
    """

    def _positive_roots_in_greater_than_zero() -> List[Any]:
        roots: List[Any] = []
        for eig, eig_roots in reduction_data.grading.grading.items():
            if int(eig) > 0:
                roots.extend(eig_roots)
        roots.sort(key=str)
        return roots

    def _simple_root_index(root: Any) -> Optional[int]:
        # 获取根的单项式系数，表示根在单根基下的分解
        # 例如：根 α = 2α₁ + 3α₂ 的单项式系数为 {1: 2, 2: 3}
        mc = root.monomial_coefficients()
        # 提取非零系数，将索引和系数转换为整数
        nonzero = [(int(i), int(c)) for i, c in mc.items() if int(c) != 0]
        # 判断是否为单根：只有一个非零系数且该系数为1
        # 例如：α₁ = 1·α₁，单项式系数为 {1: 1}
        if len(nonzero) == 1 and nonzero[0][1] == 1:
            # 返回单根的索引（在Dynkin图中的位置）
            return nonzero[0][0]
        # 如果不是单根，返回None
        return None

    def _chi_simple_support() -> Dict[int, int]:
        # f lives in negative root spaces. For each root in support of f,
        # map to the corresponding positive root and keep simple ones as
        # character support chi_i = 1.
        support: Dict[int, int] = {}
        f_mc = reduction_data.triple.f.monomial_coefficients()
        for neg_root in f_mc.keys():
            pos_root = -neg_root
            idx = _simple_root_index(pos_root)
            if idx is not None:
                support[idx] = 1
        return support

    def _ghost_grade_multiplicities(num_root_ghosts: int, max_grade: int) -> List[int]:
        """Return truncated ghost Fock multiplicities by conformal grade.

        Phase-2a uses a mode-counting model inspired by the ghost sector in
        FKW Eq. (3.1.7): for each positive root there are two fermionic ghost
        species with mode energies ``m >= 1``. At fixed energy ``m``, each
        species can be occupied at most once, so the generating factor is

            (1 + q^m) ** (2 * num_root_ghosts).

        This function extracts coefficients up to ``q^max_grade``.
        """
        if max_grade < 0:
            raise ValueError("max_grade must be non-negative")
        coeff = [0] * (max_grade + 1)
        coeff[0] = 1

        species_count = 2 * num_root_ghosts
        for mode_energy in range(1, max_grade + 1):
            for _ in range(species_count):
                for grade in range(max_grade, mode_energy - 1, -1):
                    coeff[grade] += coeff[grade - mode_energy]
        return coeff

    def _repeat_block_diagonal(base: Matrix, times: int) -> Matrix:
        """Construct block-diagonal matrix with ``times`` copies of ``base``."""
        if times < 0:
            raise ValueError("times must be non-negative")
        if times == 0:
            return Matrix(QQ, 0, 0, [])
        if times == 1:
            return base

        r = base.nrows()
        c = base.ncols()
        out = Matrix(QQ, r * times, c * times, 0)
        for blk in range(times):
            ro = blk * r
            co = blk * c
            for i in range(r):
                for j in range(c):
                    val = base[i, j]
                    if val != 0:
                        out[ro + i, co + j] = val
        return out

    def _build_d1_base(roots: Sequence[Any]) -> Matrix:
        """Build N=0 base differential d1: C^1 -> C^2 from Lie brackets."""
        dim_c1_base = len(roots)
        pairs: List[Tuple[int, int]] = []
        for i in range(dim_c1_base):
            for j in range(i + 1, dim_c1_base):
                pairs.append((i, j))
        pair_to_row = {pair: row for row, pair in enumerate(pairs)}
        dim_c2_base = len(pairs)

        if dim_c2_base == 0:
            return Matrix(QQ, 0, dim_c1_base, [])

        d1_base = Matrix(QQ, dim_c2_base, dim_c1_base, [0] * (dim_c2_base * dim_c1_base))
        lie_alg = LieAlgebra(QQ, cartan_type=["A", rank])
        basis = lie_alg.basis()

        # FKW §2.2: structure-constant (cubic ghost) term in d_st.
        for alpha_col, alpha in enumerate(roots):
            for i in range(dim_c1_base):
                for j in range(i + 1, dim_c1_base):
                    beta = roots[i]
                    gamma = roots[j]
                    bracket = lie_alg.bracket(basis[beta], basis[gamma])
                    coeff = bracket.monomial_coefficients().get(alpha, 0)
                    if coeff != 0:
                        row = pair_to_row[(i, j)]
                        d1_base[row, alpha_col] -= coeff
        return d1_base

    # Phase 2a: compute BRST cohomology for all blocks 0 <= N <= Nmax in one pass.
    def _compute_truncated_brs_cohomology(
        labels: Sequence[int],
        max_grade: int,
        roots: Sequence[Any],
        d1_base: Matrix,
        rank_d1_base: int,
    ) -> Tuple[int, int, Tuple[Tuple[int, int], ...]]:
        if max_grade < 0:
            raise ValueError("cohomology_max_grade must be non-negative")

        # Base N=0 BRST block dimensions.
        dim_c0_base = 1
        dim_c1_base = len(roots)

        if dim_c1_base == 0:
            by_grade = tuple((1 if n == 0 else 0, 0) for n in range(max_grade + 1))
            return sum(v[0] for v in by_grade), 0, by_grade

        chi_support = _chi_simple_support()

        d0_entries: List[Fraction] = []
        for root in roots:
            idx = _simple_root_index(root)
            if idx is None:
                d0_entries.append(Fraction(0, 1))
                continue
            label_val = int(labels[idx - 1]) if 1 <= idx <= rank else 0
            coeff = Fraction(label_val - chi_support.get(idx, 0), 1)
            d0_entries.append(coeff)

        d0_base = Matrix(QQ, dim_c1_base, dim_c0_base, d0_entries)

        # Base nilpotency check.
        comp0 = d1_base * d0_base
        if any(comp0[r, c] != 0 for r in range(comp0.nrows()) for c in range(comp0.ncols())):
            raise ValueError("Phase-2a base BRST nilpotency check failed: d1*d0 != 0 at N=0")

        rank_d0_base = int(d0_base.rank())
        dim_h0_base = dim_c0_base - rank_d0_base
        dim_h1_base = (dim_c1_base - rank_d1_base) - rank_d0_base
        if dim_h1_base < 0:
            raise ValueError("Invalid base BRST block: negative H^1 at N=0")

        total_h0 = 0
        total_h1 = 0
        by_grade_list: List[Tuple[int, int]] = []
        ghost_mult = _ghost_grade_multiplicities(num_root_ghosts=dim_c1_base, max_grade=max_grade)

        for grade_n in range(max_grade + 1):
            # FKW §3.2 decomposition idea: compute each L0 block separately.
            # Multiplicity comes from truncated ghost Fock counting.
            mult = ghost_mult[grade_n]

            if mult == 0:
                by_grade_list.append((0, 0))
                continue

            # Phase-2a explicit block assembly:
            #   C^0_N -> C^1_N -> C^2_N
            # by repeating the N=0 BRST block according to the grade multiplicity.
            d0_n = _repeat_block_diagonal(d0_base, mult)
            d1_n = _repeat_block_diagonal(d1_base, mult)

            # Per-grade nilpotency check, not only N=0 base check.
            comp_n = d1_n * d0_n
            if any(comp_n[r, c] != 0 for r in range(comp_n.nrows()) for c in range(comp_n.ncols())):
                raise ValueError(
                    f"Phase-2a BRST nilpotency check failed: d1*d0 != 0 at grade N={grade_n}"
                )

            dim_c0_n = dim_c0_base * mult
            dim_c1_n = dim_c1_base * mult
            rank_d0_n = int(d0_n.rank())
            rank_d1_n = int(d1_n.rank())

            # Decomposition checks: block-diagonal repetition should scale linearly.
            if rank_d0_n != mult * rank_d0_base:
                raise ValueError(
                    f"Phase-2 decomposition check failed at grade N={grade_n}: "
                    "rank(d0_N) is inconsistent with repeated base block"
                )
            if rank_d1_n != mult * rank_d1_base:
                raise ValueError(
                    f"Phase-2 decomposition check failed at grade N={grade_n}: "
                    "rank(d1_N) is inconsistent with repeated base block"
                )

            h0_n = dim_c0_n - rank_d0_n
            h1_n = (dim_c1_n - rank_d1_n) - rank_d0_n
            if h1_n < 0:
                raise ValueError(f"Invalid BRST block at grade N={grade_n}: negative H^1")

            expected_h0_n = mult * dim_h0_base
            expected_h1_n = mult * dim_h1_base
            if h0_n != expected_h0_n or h1_n != expected_h1_n:
                raise ValueError(
                    f"Phase-2 decomposition check failed at grade N={grade_n}: "
                    "computed H^q_N does not match repeated base block"
                )

            by_grade_list.append((h0_n, h1_n))
            total_h0 += h0_n
            total_h1 += h1_n

        return total_h0, total_h1, tuple(by_grade_list)

    roots = _positive_roots_in_greater_than_zero()
    d1_base = _build_d1_base(roots)
    rank_d1_base = int(d1_base.rank())

    candidates: List[AdmissibleWeightCandidate] = []
    for labels in _dominant_finite_dynkin_labels(rank=rank, affine_level=affine_level):
        total = sum(labels)
        lambda0 = affine_level - total
        grade = sum(labels[i] * int(weighted_dynkin[i]) for i in range(rank))
        h0_dim, h1_dim, by_grade = _compute_truncated_brs_cohomology(
            labels,
            cohomology_max_grade,
            roots,
            d1_base,
            rank_d1_base,
        )

        # Keep non-vacuum reduced modules with nontrivial N=0 BRST cohomology.
        n0_h0, n0_h1 = by_grade[0] if by_grade else (0, 0)
        survives = total > 0 and (n0_h0 + n0_h1) > 0

        if survives and nondegeneracy_conditions:
            affine_labels = (Fraction(lambda0, 1),) + tuple(Fraction(v, 1) for v in labels)
            survives = is_nondegenerate_affine_dynkin(affine_labels, nondegeneracy_conditions)

        candidates.append(
            AdmissibleWeightCandidate(
                finite_dynkin_labels=labels,
                lambda0_coeff=lambda0,
                affine_level=affine_level,
                reduction_grade=grade,
                cohomology_computed_upto_grade=cohomology_max_grade,
                cohomology_h0_dim=h0_dim,
                cohomology_h1_dim=h1_dim,
                cohomology_by_grade=by_grade,
                cohomology_hq_by_grade=tuple((n, by_grade[n]) for n in range(len(by_grade))),
                survives_reduction=survives,
            )
        )
    return candidates


def _compute_reduced_conformal_weight(
    rank: int,
    labels: Sequence[int],
    level: FractionalLevel,
    reduction_data: ReductionData,
) -> Fraction:
    """Compute reduced highest-weight conformal weight (MVP explicit policy)."""
    k_plus_h = _to_fraction(level.k_plus_h_vee)
    casimir_part = _finite_weight_norm_A(labels=labels, rank=rank) / (2 * k_plus_h)

    lambda_h = Fraction(0, 1)
    for i, label in enumerate(labels):
        lambda_h += Fraction(int(label), 1) * reduction_data.h_coroot_coeffs[i]

    return casimir_part - lambda_h / 2


def quantum_ds_reduction_workflow(
    cartan_type: str,
    rank: int,
    partition: Sequence[int],
    level_input: Union[PrincipalAdmissibleLevelInput, FractionalLevel],
    cohomology_max_grade: int = 0,
    nondegeneracy_conditions: Optional[Sequence[LinearNonDegeneracyCondition]] = None,
) -> QuantumDSReductionResult:
    """Run MVP quantum DS reduction workflow for partition-labelled ``sl_n`` input."""
    if not isinstance(cohomology_max_grade, int):
        raise ValueError("cohomology_max_grade must be an integer")

    letter = str(cartan_type).upper()
    if letter != "A":
        raise NotImplementedError("MVP DS workflow currently supports type A only")

    n = int(rank) + 1
    if sum(int(x) for x in partition) != n:
        raise ValueError(f"For A_{rank}, partition must sum to {n}")

    if isinstance(level_input, FractionalLevel):
        level = level_input
        level_cartan = tuple(level.cartan_type)
        if len(level_cartan) != 3:
            raise ValueError(
                "FractionalLevel cartan_type must be affine [letter, rank, 1] for MVP DS workflow"
            )

        level_letter = str(level_cartan[0]).upper()
        level_rank = int(level_cartan[1])
        level_affine_mark = int(level_cartan[2])
        if (level_letter, level_rank, level_affine_mark) != (letter, int(rank), 1):
            raise ValueError(
                "FractionalLevel cartan_type does not match requested workflow input "
                f"({letter}, {int(rank)}, 1)"
            )
    else:
        level = FractionalLevel([letter, rank, 1], p=level_input.p, u=level_input.u)

    reduction_data = _build_reduction_data(letter, int(rank), partition)

    affine_level = int(level.p) - int(level.h_vee)
    candidates = _compute_reduction_candidates(
        rank=int(rank),
        affine_level=affine_level,
        weighted_dynkin=reduction_data.weighted_dynkin,
        reduction_data=reduction_data,
        cohomology_max_grade=cohomology_max_grade,
        nondegeneracy_conditions=nondegeneracy_conditions,
    )

    k = _to_fraction(level.level)
    central_charge = _principal_central_charge_sl_n(n=n, k=k)

    modules: List[ReducedModuleData] = []
    for cand in candidates:
        h_red = _compute_reduced_conformal_weight(
            rank=int(rank),
            labels=cand.finite_dynkin_labels,
            level=level,
            reduction_data=reduction_data,
        )
        modules.append(
            ReducedModuleData(
                candidate=cand,
                central_charge=central_charge,
                conformal_weight=h_red,
            )
        )

    surviving_candidates = [c for c in candidates if c.survives_reduction]
    surviving_modules = [m for m in modules if m.candidate.survives_reduction]

    return QuantumDSReductionResult(
        cartan_type=letter,
        rank=int(rank),
        partition=tuple(int(x) for x in partition),
        level=level,
        reduction_data=reduction_data,
        all_candidates=candidates,
        surviving_candidates=surviving_candidates,
        all_modules=modules,
        surviving_modules=surviving_modules,
    )


__all__ = [
    "PrincipalAdmissibleLevelInput",
    "ReductionData",
    "LinearNonDegeneracyCondition",
    "is_nondegenerate_affine_dynkin",
    "AdmissibleWeightCandidate",
    "ReducedModuleData",
    "CharacterEvaluation",
    "QuantumDSReductionResult",
    "quantum_ds_reduction_workflow",
]
