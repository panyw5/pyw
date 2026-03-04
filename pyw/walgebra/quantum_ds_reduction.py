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

from sage.all import CartanType, LieAlgebra, Matrix, QQ

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
    cohomology_h0_dim: int
    cohomology_h1_dim: int
    survives_reduction: bool


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

    def _compute_n0_brs_cohomology(labels: Sequence[int]) -> Tuple[int, int]:
        # Phase 1.5 keeps only the N=0 block of the BRST complex from FKW §2.3:
        #   C(M) = M \otimes \Lambda, with differential d = d_st + p
        # Here we model the first three ghost-number pieces C^0 -> C^1 -> C^2.
        roots = _positive_roots_in_greater_than_zero()
        dim_c0 = 1
        dim_c1 = len(roots)

        if dim_c1 == 0:
            return 1, 0

        # C^2 basis: ordered pairs (i, j), i < j
        pairs: List[Tuple[int, int]] = []
        for i in range(dim_c1):
            for j in range(i + 1, dim_c1):
                pairs.append((i, j))
        pair_to_row = {pair: row for row, pair in enumerate(pairs)}
        dim_c2 = len(pairs)

        chi_support = _chi_simple_support()

        # d0: C^0 -> C^1.
        # The shift by chi corresponds to the character term p_\pm in FKW
        # (see Eq. (2.1.1) and d_\pm = d_st^\pm + p_\pm around Eq. (2.2.1)- (2.2.2) block).
        d0_entries: List[Fraction] = []
        for root in roots:
            idx = _simple_root_index(root)
            if idx is None:
                d0_entries.append(Fraction(0, 1))
                continue

            # Type A uses simple-root indices 1..rank in Dynkin labels tuple.
            label_val = int(labels[idx - 1]) if 1 <= idx <= rank else 0
            coeff = Fraction(label_val - chi_support.get(idx, 0), 1)
            d0_entries.append(coeff)

        d0 = Matrix(QQ, dim_c1, dim_c0, d0_entries)

        # d1: C^1 -> C^2.
        # We now use actual Lie brackets [e_beta, e_gamma] = sum c_{beta,gamma}^alpha e_alpha
        # from Sage's Chevalley basis, matching the structure-constant term in
        # d_st (FKW Eq. for d_st^\pm(z), cubic ghost term).
        if dim_c2 == 0:
            d1 = Matrix(QQ, 0, dim_c1, [])
        else:
            d1 = Matrix(QQ, dim_c2, dim_c1, [0] * (dim_c2 * dim_c1))
            lie_alg = LieAlgebra(QQ, cartan_type=["A", rank])
            basis = lie_alg.basis()

            for alpha_col, alpha in enumerate(roots):
                for i in range(dim_c1):
                    for j in range(i + 1, dim_c1):
                        beta = roots[i]
                        gamma = roots[j]
                        bracket = lie_alg.bracket(basis[beta], basis[gamma])
                        coeff = bracket.monomial_coefficients().get(alpha, 0)
                        if coeff != 0:
                            row = pair_to_row[(i, j)]
                            d1[row, alpha_col] -= coeff

        # Nilpotency check for this finite block: d1 * d0 must vanish.
        # This is the matrix-level shadow of d^2 = 0 in BRST theory.
        composed = d1 * d0
        if any(
            composed[r, c] != 0 for r in range(composed.nrows()) for c in range(composed.ncols())
        ):
            raise ValueError("Phase-1.5 BRST nilpotency check failed: d1*d0 != 0")

        rank_d0 = int(d0.rank())
        rank_d1 = int(d1.rank())

        dim_ker_d0 = dim_c0 - rank_d0
        dim_ker_d1 = dim_c1 - rank_d1
        dim_h0 = dim_ker_d0
        dim_h1 = dim_ker_d1 - rank_d0

        if dim_h1 < 0:
            raise ValueError("Invalid BRST complex block: im(d0) not contained in ker(d1)")

        return int(dim_h0), int(dim_h1)

    candidates: List[AdmissibleWeightCandidate] = []
    for labels in _dominant_finite_dynkin_labels(rank=rank, affine_level=affine_level):
        total = sum(labels)
        lambda0 = affine_level - total
        grade = sum(labels[i] * int(weighted_dynkin[i]) for i in range(rank))
        h0_dim, h1_dim = _compute_n0_brs_cohomology(labels)

        # Keep non-vacuum reduced modules with nontrivial N=0 BRST cohomology.
        survives = total > 0 and (h0_dim + h1_dim) > 0

        if survives and nondegeneracy_conditions:
            affine_labels = (Fraction(lambda0, 1),) + tuple(Fraction(v, 1) for v in labels)
            survives = is_nondegenerate_affine_dynkin(affine_labels, nondegeneracy_conditions)

        candidates.append(
            AdmissibleWeightCandidate(
                finite_dynkin_labels=labels,
                lambda0_coeff=lambda0,
                affine_level=affine_level,
                reduction_grade=grade,
                cohomology_h0_dim=h0_dim,
                cohomology_h1_dim=h1_dim,
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
    nondegeneracy_conditions: Optional[Sequence[LinearNonDegeneracyCondition]] = None,
) -> QuantumDSReductionResult:
    """Run MVP quantum DS reduction workflow for partition-labelled ``sl_n`` input."""
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
