from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from itertools import product
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Tuple

from sage.all import IntegrableRepresentation as SageIntegrableRepresentation
from sage.all import Integer, QQ, SR, ZZ, binomial, matrix, prod, var, vector

if TYPE_CHECKING:
    from .affine_lie_algebra import AffineLieAlgebra
    from .affine_weight import AffineWeight


def _element_word_list(element: Any) -> list[int]:
    if hasattr(element, "reduced_word_list"):
        return list(element.reduced_word_list())
    return [int(i) for i in element.reduced_word()]


def _to_extended_affine_highest_weight(algebra: "AffineLieAlgebra", highest_weight: Any) -> Any:
    """Rebuild an affine highest weight in Sage's extended affine lattice."""
    if not algebra.is_affine:
        return highest_weight

    extended_affine_weight_lattice = algebra.affine_weight_lattice_sage()
    extended_affine_fundamental_weights = extended_affine_weight_lattice.fundamental_weights()
    rebuilt_highest_weight = extended_affine_weight_lattice.zero()
    for index, coefficient in highest_weight.monomial_coefficients().items():
        integer_coefficient = Integer(coefficient)
        if integer_coefficient != 0:
            rebuilt_highest_weight += integer_coefficient * extended_affine_fundamental_weights[index]
    return rebuilt_highest_weight


def _apply_affine_element_to_weight(
    algebra: "AffineLieAlgebra", element: Any, weight: "AffineWeight"
) -> "AffineWeight":
    from .affine_weight import AffineWeight

    try:
        domain_weight = weight.to_sagemath()
        acted_weight = element.action(domain_weight)
        acted_vector = acted_weight.to_vector()
        acted_grade = QQ(acted_vector[-1]) if len(acted_vector) > 0 else QQ(0)
        return AffineWeight.from_sagemath(algebra, acted_weight, grade=acted_grade)
    except Exception:
        pass

    semidirect = algebra.affine_weyl_group()
    word = tuple(_element_word_list(element))
    return semidirect.from_word(word).action(weight)


def _finite_coroot_gram_matrix(algebra: "AffineLieAlgebra", idxs: List[int]) -> List[List[Any]]:
    semidirect = algebra.affine_weyl_group()
    basis = semidirect._finite_coroot_space.simple_roots()
    ambient = {i: basis[i].to_ambient() for i in idxs}
    return [[QQ(ambient[i].inner_product(ambient[j])) for j in idxs] for i in idxs]


def _ceil_sqrt_qq(value: Any) -> int:
    target = QQ(value)
    if target <= 0:
        return 0

    lower = 0
    upper = 1
    while QQ(upper) * QQ(upper) < target:
        lower = upper
        upper *= 2

    while lower + 1 < upper:
        mid = (lower + upper) // 2
        if QQ(mid) * QQ(mid) >= target:
            upper = mid
        else:
            lower = mid

    return upper


def _translation_coefficient_radius(
    *,
    level: Any,
    linear_coeffs: List[Any],
    gram: List[List[Any]],
    max_neg_shift: Any,
) -> int:
    k = QQ(level)
    b = QQ(_ceil_sqrt_qq(sum(QQ(d) * QQ(d) for d in linear_coeffs)))
    G = matrix(QQ, gram)
    G_inv = G.inverse()
    frob_sq = sum(QQ(v) * QQ(v) for v in G_inv.list())
    lam_lower = QQ(1) / QQ(_ceil_sqrt_qq(frob_sq))

    if lam_lower <= 0:
        raise ValueError("Failed to obtain a positive coercive bound for translation enumeration")

    abs_k = abs(k)
    a = abs_k * lam_lower / QQ(2)
    if a <= 0:
        raise ValueError("Failed to derive a positive quadratic bound for translation enumeration")

    if k > 0:
        disc = b * b + QQ(4) * a * QQ(max_neg_shift)
        r_real = (b + QQ(_ceil_sqrt_qq(disc))) / (QQ(2) * a)
        return max(0, int(r_real) + 1)

    r_real = b / a
    return max(0, int(r_real) + 1)


def _minus_delta_n(
    *,
    level: Any,
    linear_coeffs: List[Any],
    gram: List[List[Any]],
    coeffs: tuple[int, ...],
) -> Any:
    m = [QQ(c) for c in coeffs]
    linear = sum(QQ(d) * mi for d, mi in zip(linear_coeffs, m))
    quadratic = QQ(0)
    for i, mi in enumerate(m):
        for j, mj in enumerate(m):
            quadratic += mi * QQ(gram[i][j]) * mj
    return linear + QQ(level) * quadratic / QQ(2)


def _translations_by_n_shift_bnb_impl(
    algebra: "AffineLieAlgebra",
    weight: "AffineWeight",
    *,
    order: int,
    max_neg_shift: Optional[Any] = None,
    return_stats: bool = False,
) -> Any:
    semidirect = algebra.affine_weyl_group()
    idxs = [int(i) for i in semidirect._finite_coroot_space.index_set()]
    basis = semidirect._finite_coroot_space.simple_roots()

    if max_neg_shift is None:
        max_neg_shift_value = QQ(order) + QQ(weight.grade)
    else:
        max_neg_shift_value = QQ(max_neg_shift)
    if max_neg_shift_value < 0:
        return ({"translations": [], "stats": {}} if return_stats else [])

    level = QQ(weight.level)
    if level == 0:
        raise ValueError("Translation enumeration by n-shift requires non-zero level")

    linear_coeffs = [QQ(weight.dynkin_labels().get(i, 0)) for i in idxs]
    gram = _finite_coroot_gram_matrix(algebra, idxs)
    radius = _translation_coefficient_radius(
        level=level,
        linear_coeffs=linear_coeffs,
        gram=gram,
        max_neg_shift=max_neg_shift_value,
    )

    n = len(idxs)
    if n == 0:
        t0 = semidirect.translation(semidirect._zero_beta)
        return ({"translations": [t0], "stats": {"radius": 0}} if return_stats else [t0])

    center: list[QQ]
    try:
        G = matrix(QQ, gram)
        d_vec = vector(QQ, linear_coeffs)
        center_vec = -(G.solve_right(d_vec)) / level
        center = [QQ(center_vec[i]) for i in range(n)]
    except Exception:
        center = [QQ(0) for _ in range(n)]

    coeffs = [0 for _ in range(n)]
    selected: dict[tuple[int, ...], Any] = {}

    stats = {
        "radius": int(radius),
        "dimension": int(n),
        "box_points": int((2 * int(radius) + 1) ** int(n)),
        "visited_leaves": 0,
        "pruned_branches": 0,
        "accepted": 0,
    }

    gram_qq = [[QQ(gram[i][j]) for j in range(n)] for i in range(n)]
    tail_quad_lower = [QQ(0) for _ in range(n + 1)]
    tail_quad_upper = [QQ(0) for _ in range(n + 1)]
    R = QQ(radius)
    for depth in range(n - 1, -1, -1):
        lower = tail_quad_lower[depth + 1]
        upper = tail_quad_upper[depth + 1]
        qii = level * gram_qq[depth][depth] / QQ(2)
        term = qii * R * R
        if term >= 0:
            upper += term
        else:
            lower += term
        for j in range(depth + 1, n):
            band = abs(level * gram_qq[depth][j]) * R * R
            lower -= band
            upper += band
        tail_quad_lower[depth] = lower
        tail_quad_upper[depth] = upper

    current_const = QQ(0)
    current_b = [QQ(v) for v in linear_coeffs]

    def _partial_bounds(depth: int) -> tuple[Any, Any]:
        if depth == n:
            return current_const, current_const

        lower = QQ(current_const) + tail_quad_lower[depth]
        upper = QQ(current_const) + tail_quad_upper[depth]

        for i in range(depth, n):
            delta = abs(current_b[i]) * R
            lower -= delta
            upper += delta

        return lower, upper

    ordered_values_by_dim: list[list[int]] = []
    base_values = list(range(-radius, radius + 1))
    for i in range(n):
        vals = list(base_values)
        c = center[i]
        vals.sort(key=lambda x: (abs(QQ(x) - c), abs(x), x))
        ordered_values_by_dim.append(vals)

    def _dfs(depth: int, prefix_norm_sq: int) -> None:
        nonlocal current_const, current_b
        if prefix_norm_sq > radius * radius:
            stats["pruned_branches"] += 1
            return

        low, high = _partial_bounds(depth)
        if high < 0 or low > max_neg_shift_value:
            stats["pruned_branches"] += 1
            return

        if depth == n:
            stats["visited_leaves"] += 1
            neg_shift = _minus_delta_n(
                level=level,
                linear_coeffs=linear_coeffs,
                gram=gram,
                coeffs=tuple(coeffs),
            )
            if QQ(0) <= neg_shift <= max_neg_shift_value:
                beta = semidirect._zero_beta
                for idx, c in zip(idxs, coeffs):
                    if c:
                        beta += int(c) * basis[idx]
                key = tuple(int(c) for c in coeffs)
                selected[key] = semidirect.translation(beta)
                stats["accepted"] += 1
            return

        for value in ordered_values_by_dim[depth]:
            coeffs[depth] = int(value)
            value_qq = QQ(value)
            old_const = current_const
            old_b = list(current_b)

            current_const = old_const + current_b[depth] * value_qq + level * gram_qq[depth][depth] * value_qq * value_qq / QQ(2)
            for j in range(depth + 1, n):
                current_b[j] = old_b[j] + level * gram_qq[j][depth] * value_qq

            _dfs(depth + 1, prefix_norm_sq + int(value) * int(value))
            current_const = old_const
            current_b = old_b
        coeffs[depth] = 0

    _dfs(0, 0)
    out = [selected[key] for key in sorted(selected.keys())]
    if return_stats:
        return {"translations": out, "stats": stats}
    return out


def _translations_by_n_shift_impl(
    algebra: "AffineLieAlgebra",
    weight: "AffineWeight",
    *,
    order: int,
    max_neg_shift: Optional[Any] = None,
) -> List[Any]:
    semidirect = algebra.affine_weyl_group()
    idxs = [int(i) for i in semidirect._finite_coroot_space.index_set()]
    basis = semidirect._finite_coroot_space.simple_roots()

    if max_neg_shift is None:
        max_neg_shift_value = QQ(order) + QQ(weight.grade)
    else:
        max_neg_shift_value = QQ(max_neg_shift)
    if max_neg_shift_value < 0:
        return []

    if QQ(weight.level) == 0:
        raise ValueError("Translation enumeration by n-shift requires non-zero level")

    linear_coeffs = [QQ(weight.dynkin_labels().get(i, 0)) for i in idxs]
    gram = _finite_coroot_gram_matrix(algebra, idxs)
    radius = _translation_coefficient_radius(
        level=QQ(weight.level),
        linear_coeffs=linear_coeffs,
        gram=gram,
        max_neg_shift=max_neg_shift_value,
    )

    dimension = len(idxs)
    box_points = (2 * int(radius) + 1) ** int(dimension)
    if box_points > 2_000_000:
        return _translations_by_n_shift_bnb_impl(
            algebra,
            weight,
            order=order,
            max_neg_shift=max_neg_shift_value,
            return_stats=False,
        )

    ranges = [range(-radius, radius + 1) for _ in idxs]
    selected: dict[tuple[int, ...], Any] = {}
    for coeffs in product(*ranges):
        neg_shift = _minus_delta_n(
            level=QQ(weight.level),
            linear_coeffs=linear_coeffs,
            gram=gram,
            coeffs=coeffs,
        )
        if neg_shift < 0 or neg_shift > max_neg_shift_value:
            continue

        beta = semidirect._zero_beta
        for i, coeff in zip(idxs, coeffs):
            if coeff:
                beta += int(coeff) * basis[i]
        selected[tuple(int(coeff) for coeff in coeffs)] = semidirect.translation(beta)

    return [selected[key] for key in sorted(selected.keys())]


def _build_W_affine_as_words_direct_impl(
    algebra: "AffineLieAlgebra",
    normalized_translation_vectors: Iterable[Any],
    *,
    finite_affine_elements_cache: Optional[List[Any]] = None,
) -> List[Any]:
    affine_weyl_group = algebra.affine_weyl_group()
    affine_weyl_group_sage = algebra.affine_weyl_group_sage()

    finite_affine_elements = finite_affine_elements_cache
    if finite_affine_elements is None:
        finite_words = [
            tuple(int(i) for i in w.reduced_word())
            for w in list(affine_weyl_group._finite_weyl_group)
        ]
        finite_affine_elements = [
            affine_weyl_group_sage.from_reduced_word(list(word))
            for word in finite_words
        ]

    candidates: List[Any] = []
    for beta in normalized_translation_vectors:
        translation_word = tuple(int(i) for i in affine_weyl_group.translation_word_list(beta))
        translation_sage = affine_weyl_group_sage.from_reduced_word(list(translation_word))
        for finite_affine in finite_affine_elements:
            candidates.append(finite_affine * translation_sage)

    return candidates


@dataclass(frozen=True)
class KLNumeratorTerm:
    representative: Any
    weight: "AffineWeight"
    coefficient: Any


@dataclass
class KazhdanLusztigData:
    algebra: Any
    lambda_hat: "AffineWeight"
    Lambda_hat: "AffineWeight"
    w_to_lambda: Any
    order: int
    translations: List[Any]
    W_affine_as_words: List[Any]
    stabilizer_candidates: List[Any]
    quotient_representatives: List[Any]

    @property
    def manual_translations(self) -> List[Any]:
        """Backward-compatible alias for older call sites/tests."""
        return self.translations

    @property
    def integral_case(self) -> bool:
        labels = self.Lambda_hat.dynkin_labels()
        return all(value in ZZ for value in labels.values())

    def apply(self, element: Any, weight: "AffineWeight") -> "AffineWeight":
        semidirect = self.algebra.affine_weyl_group()
        word = tuple(_element_word_list(element))
        return semidirect.from_word(word).action(weight)

    def quotient_weight(self, representative: Any) -> "AffineWeight":
        rho_hat = self.algebra.affine_rho()
        return self.apply(representative, self.Lambda_hat + rho_hat) - rho_hat

    def weyl_to_be_summed(self) -> List[Any]:
        from .bruhat import BruhatOrder

        if not self.W_affine_as_words:
            return []
        bruhat = BruhatOrder(self.W_affine_as_words[0].parent())
        lower = self.w_to_lambda
        rep_map = {
            tuple(_element_word_list(w)): w for w in self.quotient_representatives
        }
        if not rep_map:
            return []

        max_length = max(int(w.length()) for w in self.quotient_representatives)
        queue = deque([lower])
        visited = {tuple(_element_word_list(lower))}
        result_map: Dict[Tuple[int, ...], Any] = {}

        while queue:
            current = queue.popleft()
            current_key = tuple(_element_word_list(current))
            current_length = int(current.length())

            if current_key in rep_map:
                result_map[current_key] = rep_map[current_key]

            if current_length >= max_length:
                continue

            for nxt in current.bruhat_upper_covers():
                nxt_key = tuple(_element_word_list(nxt))
                if nxt_key not in visited:
                    visited.add(nxt_key)
                    queue.append(nxt)

        return sorted(result_map.values(), key=lambda w: (int(w.length()), tuple(_element_word_list(w))))


@dataclass
class FormalCharacter:
    coefficients: Dict[int, Any] = field(default_factory=dict)
    max_grade: int = 10
    algebra: Optional["AffineLieAlgebra"] = None

    def __post_init__(self) -> None:
        self.coefficients = {
            grade: coeff
            for grade, coeff in self.coefficients.items()
            if coeff != 0 and grade <= self.max_grade
        }

    def __getitem__(self, grade: int) -> Any:
        return self.coefficients.get(grade, 0)

    def __setitem__(self, grade: int, value: Any) -> None:
        if value == 0:
            self.coefficients.pop(grade, None)
            return
        if grade <= self.max_grade:
            self.coefficients[grade] = value

    def __iter__(self) -> Iterator[Tuple[int, Any]]:
        for grade in sorted(self.coefficients):
            yield grade, self.coefficients[grade]

    def __len__(self) -> int:
        return len(self.coefficients)

    def __add__(self, other: "FormalCharacter") -> "FormalCharacter":
        if not isinstance(other, FormalCharacter):
            return NotImplemented

        max_grade = max(self.max_grade, other.max_grade)
        result: Dict[int, Any] = {}
        for grade in set(self.coefficients) | set(other.coefficients):
            if grade > max_grade:
                continue
            coeff = self[grade] + other[grade]
            if coeff != 0:
                result[grade] = coeff
        return FormalCharacter(result, max_grade=max_grade, algebra=self.algebra or other.algebra)

    def __sub__(self, other: "FormalCharacter") -> "FormalCharacter":
        if not isinstance(other, FormalCharacter):
            return NotImplemented

        max_grade = max(self.max_grade, other.max_grade)
        result: Dict[int, Any] = {}
        for grade in set(self.coefficients) | set(other.coefficients):
            if grade > max_grade:
                continue
            coeff = self[grade] - other[grade]
            if coeff != 0:
                result[grade] = coeff
        return FormalCharacter(result, max_grade=max_grade, algebra=self.algebra or other.algebra)

    def __mul__(self, other: Any) -> "FormalCharacter":
        if isinstance(other, FormalCharacter):
            max_grade = min(self.max_grade, other.max_grade)
            result: Dict[int, Any] = {}
            for grade_1, coeff_1 in self.coefficients.items():
                for grade_2, coeff_2 in other.coefficients.items():
                    grade = grade_1 + grade_2
                    if grade > max_grade:
                        continue
                    result[grade] = result.get(grade, 0) + coeff_1 * coeff_2
            return FormalCharacter(
                result, max_grade=max_grade, algebra=self.algebra or other.algebra
            )

        if isinstance(other, (int, Integer)) or hasattr(other, "__rmul__"):
            result = {
                grade: coeff * other
                for grade, coeff in self.coefficients.items()
                if coeff * other != 0
            }
            return FormalCharacter(result, max_grade=self.max_grade, algebra=self.algebra)

        return NotImplemented

    def __rmul__(self, scalar: Any) -> "FormalCharacter":
        return self.__mul__(scalar)

    def __neg__(self) -> "FormalCharacter":
        return self * (-1)

    @property
    def leading_grade(self) -> Optional[int]:
        if not self.coefficients:
            return None
        return min(self.coefficients)

    @property
    def leading_coefficient(self) -> Any:
        leading_grade = self.leading_grade
        if leading_grade is None:
            return 0
        return self.coefficients[leading_grade]

    def truncate(self, new_max: int) -> "FormalCharacter":
        return FormalCharacter(
            {grade: coeff for grade, coeff in self.coefficients.items() if grade <= new_max},
            max_grade=new_max,
            algebra=self.algebra,
        )

    def shift(self, delta: int) -> "FormalCharacter":
        shifted = {grade + delta: coeff for grade, coeff in self.coefficients.items()}
        return FormalCharacter(shifted, max_grade=self.max_grade + delta, algebra=self.algebra)

    def to_dict(self) -> Dict[int, Any]:
        return dict(self.coefficients)

    def to_list(self, up_to: Optional[int] = None) -> list[Any]:
        if up_to is None:
            up_to = self.max_grade
        return [self[grade] for grade in range(up_to + 1)]


class WeylKacDenominator:
    def __init__(self, algebra: "AffineLieAlgebra") -> None:
        self.algebra = algebra
        self._inverse_cache: Dict[int, FormalCharacter] = {}

    def inverse(self, max_grade: int = 10) -> FormalCharacter:
        cached = self._inverse_cache.get(max_grade)
        if cached is not None:
            return cached

        inverse = self._compute_inverse_product(max_grade)
        self._inverse_cache[max_grade] = inverse
        return inverse

    def _compute_inverse_product(self, max_grade: int) -> FormalCharacter:
        rank = int(getattr(self.algebra, "rank", 1) or 1)
        coeffs: Dict[int, Any] = {0: QQ(1)}

        for step in range(1, max_grade + 1):
            updated: Dict[int, Any] = {}
            for base_grade, base_coeff in coeffs.items():
                max_power = (max_grade - base_grade) // step
                for power in range(max_power + 1):
                    grade = base_grade + power * step
                    contribution = base_coeff * binomial(rank + power - 1, power)
                    updated[grade] = updated.get(grade, 0) + contribution
            coeffs = updated

        return FormalCharacter(coeffs, max_grade=max_grade, algebra=self.algebra)


class VermaCharacter:
    def __init__(self, algebra: "AffineLieAlgebra", weight: "AffineWeight") -> None:
        self.algebra = algebra
        self.weight = weight
        self._denominator = WeylKacDenominator(algebra)

    def character(self, max_grade: int = 10) -> FormalCharacter:
        weight_grade = int(getattr(self.weight, "grade", 0))
        inverse = self._denominator.inverse(max_grade + abs(weight_grade))
        return inverse.shift(-weight_grade).truncate(max_grade)


class IntegrableModuleCharacter:
    """Integrable affine highest-weight module wrapper with legacy character assembly.

    This wraps Sage's ``IntegrableRepresentation`` and also exposes the legacy
    ``CharacterOfIntegrableModule`` algorithm from ``demos/Algebra.py`` as the
    instance method ``character(order)``.
    """

    def __init__(self, highest_weight: Any) -> None:
        from .affine_lie_algebra import AffineLieAlgebra
        from .affine_weight import AffineWeight

        if isinstance(highest_weight, AffineWeight):
            self.algebra = highest_weight.algebra
            self._highest_weight_affine = highest_weight
            sage_weight = _to_extended_affine_highest_weight(
                self.algebra,
                highest_weight.to_sagemath(),
            )
        else:
            cartan_type = list(highest_weight.parent().cartan_type())
            self.algebra = AffineLieAlgebra(cartan_type)
            sage_weight = _to_extended_affine_highest_weight(self.algebra, highest_weight)
        self._highest_weight_affine = AffineWeight.from_sagemath(self.algebra, sage_weight, grade=0)

        self._highest_weight = highest_weight
        self._sage_representation = SageIntegrableRepresentation(sage_weight)

    @property
    def sage(self) -> Any:
        return self._sage_representation

    @property
    def highest_weight_affine(self) -> "AffineWeight":
        return self._highest_weight_affine

    def highest_weight(self) -> Any:
        return self._sage_representation.highest_weight()

    def dominant_maximal_weights(self) -> list[Any]:
        return list(self._sage_representation.dominant_maximal_weights())

    def strings(self, depth: int = 12) -> dict[Any, list[Any]]:
        d = int(depth)
        if d <= 0:
            raise ValueError("depth must be a positive integer")
        raw = self._sage_representation.strings(d)
        return {weight: list(values) for weight, values in raw.items()}

    def multiplicity(self, index_tuple: tuple[int, ...]) -> Any:
        return self._sage_representation.m(index_tuple)

    def to_weight(self, index_tuple: tuple[int, ...]) -> Any:
        return self._sage_representation.to_weight(index_tuple)

    def from_weight(self, weight: Any) -> tuple[int, ...]:
        return tuple(self._sage_representation.from_weight(weight))

    def _sage_weight_to_affine_weight(self, weight: Any) -> "AffineWeight":
        from .affine_weight import AffineWeight

        # In Sage's integrable highest-weight module coordinates, only alpha_0
        # contributes a delta-shift, so the grade is minus the alpha_0 count.
        root_coordinates = self.from_weight(weight)
        grade = -QQ(root_coordinates[0]) if root_coordinates else QQ(0)
        return AffineWeight.from_sagemath(self.algebra, weight, grade=grade)

    def _stable_strings_prefix(
        self,
        required_n_by_weight: dict[Any, int],
        *,
        max_rounds: int = 8,
    ) -> dict[Any, list[Any]]:
        required = {w: int(n) for w, n in required_n_by_weight.items() if int(n) >= 0}
        if not required:
            return {}

        depth = max(required.values()) + 1
        previous_prefix: Optional[dict[Any, tuple[Any, ...]]] = None
        last_data: Optional[dict[Any, list[Any]]] = None

        for _ in range(max_rounds):
            data = self.strings(depth)
            prefixes: dict[Any, tuple[Any, ...]] = {}
            complete = True
            for weight, nmax in required.items():
                if weight not in data:
                    complete = False
                    break
                seq = data[weight]
                if len(seq) < nmax + 1:
                    complete = False
                    break
                prefixes[weight] = tuple(seq[0 : nmax + 1])

            if complete and previous_prefix is not None and prefixes == previous_prefix:
                return data

            if complete:
                previous_prefix = prefixes
                last_data = data

            depth *= 2

        if last_data is not None:
            return last_data
        return self.strings(depth)

    def _auto_translations(self, dominant_weights: list[Any], *, order: int) -> list[Any]:
        affine_weyl_group = self.algebra.affine_weyl_group()
        translation_vectors: dict[Tuple[int, ...], Any] = {}

        for weight in dominant_weights:
            affine_weight = self._sage_weight_to_affine_weight(weight)
            for translation in _translations_by_n_shift_impl(
                self.algebra,
                affine_weight,
                order=order,
                max_neg_shift=QQ(order),
            ):
                beta = translation.translation_vector
                key = tuple(
                    int(QQ(beta.monomial_coefficients().get(i, 0)))
                    for i in affine_weyl_group._finite_coroot_space.index_set()
                )
                translation_vectors[key] = beta

        return [affine_weyl_group.translation(beta) for _, beta in sorted(translation_vectors.items())]


    def _weyl_candidates(self, translations: Iterable[Any]) -> list[Any]:
        affine_weyl_group = self.algebra.affine_weyl_group()
        vectors = affine_weyl_group._translations_to_coroots(translations=translations)
        return _build_W_affine_as_words_direct_impl(self.algebra, vectors)

    def _orbit_representatives_for_weight(
        self,
        weight: Any,
        *,
        candidates: Iterable[Any],
        order: int,
    ) -> list[tuple[Any, int]]:
        base_weight = self._sage_weight_to_affine_weight(weight)
        by_image: Dict[Tuple[Tuple[int, Any], ...], tuple[Any, "AffineWeight"]] = {}

        for w in candidates:
            acted = _apply_affine_element_to_weight(self.algebra, w, base_weight)
            key = tuple(sorted(acted.dynkin_labels().items())) + ((-1, acted.grade),)
            current = by_image.get(key)
            if current is None or int(w.length()) < int(current[0].length()):
                by_image[key] = (w, acted)

        selected: list[tuple[Any, int]] = []
        for wrep, acted in by_image.values():
            d0 = int(-QQ(acted.grade))
            if d0 > int(order):
                continue
            nmax = int(order) - d0
            if nmax >= 0:
                selected.append((wrep, nmax))

        return sorted(selected, key=lambda item: (int(item[0].length()), tuple(_element_word_list(item[0]))))

    def _character_contribution(
        self,
        weight: Any,
        representative: Any,
        string: list[Any],
        *,
        nmax: int,
        q_var: Any,
        z_vars: dict[int, Any],
    ) -> Any:
        from .affine_weight import AffineWeight

        base_weight = self._sage_weight_to_affine_weight(weight)
        delta = AffineWeight.delta(self.algebra)
        simple_roots = {
            i: AffineWeight.affine_simple_root(self.algebra, i) for i in range(1, self.algebra.rank + 1)
        }

        upper = min(int(nmax), len(string) - 1)
        if upper < 0:
            return 0

        return sum(
            [
                string[n]
                * prod(
                    [
                        z_vars[i]
                        ** _apply_affine_element_to_weight(
                            self.algebra,
                            representative,
                            base_weight - n * delta,
                        ).scalar_product(simple_roots[i])
                        for i in range(1, self.algebra.rank + 1)
                    ]
                )
                * q_var
                ** (
                    -_apply_affine_element_to_weight(
                        self.algebra,
                        representative,
                        base_weight - n * delta,
                    ).grade
                )
                for n in range(0, upper + 1)
            ]
        )

    def character(self, order: int, *, manual_translations: Optional[Iterable[Any]] = None) -> Any:
        target_order = int(order)
        if target_order < 0:
            return 0

        dominant_weights = self.dominant_maximal_weights()
        if manual_translations is None:
            manual_translations = self._auto_translations(dominant_weights, order=target_order)

        candidates = self._weyl_candidates(manual_translations)
        reps_by_weight: dict[Any, list[tuple[Any, int]]] = {}
        max_needed_depth_by_weight: dict[Any, int] = {}

        for weight in dominant_weights:
            reps = self._orbit_representatives_for_weight(weight, candidates=candidates, order=target_order)
            reps_by_weight[weight] = reps
            max_needed_depth_by_weight[weight] = max([nmax for _, nmax in reps], default=-1)

        if max(max_needed_depth_by_weight.values(), default=-1) < 0:
            return 0

        strings_data = self._stable_strings_prefix(max_needed_depth_by_weight)
        q_var = var("q")
        z_vars = {i: var(f"z{i}") for i in range(1, self.algebra.rank + 1)}

        total = SR(0)
        for weight in dominant_weights:
            if max_needed_depth_by_weight[weight] < 0:
                continue
            string = strings_data[weight]
            for representative, nmax in reps_by_weight[weight]:
                total += self._character_contribution(
                    weight,
                    representative,
                    string,
                    nmax=nmax,
                    q_var=q_var,
                    z_vars=z_vars,
                )

        return sum(total.coefficient(q_var, n) * q_var**n for n in range(0, target_order + 1))

    def __repr__(self) -> str:
        return repr(self._sage_representation)


class KazhdanLusztigCharacter:
    def __init__(self, algebra: "AffineLieAlgebra") -> None:
        from .kazhdan_lusztig import KazhdanLusztigPolynomials

        self.algebra = algebra
        self.kl = KazhdanLusztigPolynomials(algebra.affine_weyl_group_sage())
        self._finite_affine_elements_cache: Optional[List[Any]] = None

    @staticmethod
    def _ceil_sqrt_qq(value: Any) -> int:
        return _ceil_sqrt_qq(value)

    def _translations_by_n_shift(
        self,
        weight: "AffineWeight",
        *,
        order: int,
        max_neg_shift: Optional[Any] = None,
    ) -> List[Any]:
        return _translations_by_n_shift_impl(
            self.algebra,
            weight,
            order=order,
            max_neg_shift=max_neg_shift,
        )

    @staticmethod
    def _minus_Delta_n(
        *,
        level: Any,
        linear_coeffs: List[Any],
        gram: List[List[Any]],
        coeffs: tuple[int, ...],
    ) -> Any:
        return _minus_delta_n(
            level=level,
            linear_coeffs=linear_coeffs,
            gram=gram,
            coeffs=coeffs,
        )

    def _finite_coroot_gram_matrix(self, idxs: List[int]) -> List[List[Any]]:
        return _finite_coroot_gram_matrix(self.algebra, idxs)

    def _translation_coefficient_radius(
        self,
        *,
        level: Any,
        linear_coeffs: List[Any],
        gram: List[List[Any]],
        max_neg_shift: Any,
    ) -> int:
        return _translation_coefficient_radius(
            level=level,
            linear_coeffs=linear_coeffs,
            gram=gram,
            max_neg_shift=max_neg_shift,
        )

    def _translations_by_n_shift_bnb(
        self,
        weight: "AffineWeight",
        *,
        order: int,
        max_neg_shift: Optional[Any] = None,
        return_stats: bool = False,
    ) -> Any:
        return _translations_by_n_shift_bnb_impl(
            self.algebra,
            weight,
            order=order,
            max_neg_shift=max_neg_shift,
            return_stats=return_stats,
        )

    def _translation_neg_shift(self, weight: "AffineWeight", beta: Any) -> Any:
        translated = self.algebra.affine_weyl_group().translation(beta).action(weight)
        return weight.grade - translated.grade

    def _find_dominant_Lambda(
        self,
        lambda_hat: "AffineWeight",
        *,
        max_steps: int = 1000,
    ) -> Tuple["AffineWeight", Any]:
        rho_hat = self.algebra.affine_rho()
        weyl_group = self.kl.weyl_group
        current = lambda_hat + rho_hat
        current_w = weyl_group.one()
        visited: set[tuple[tuple[int, Any], ...]] = set()

        for _ in range(max_steps):
            labels = current.dynkin_labels()
            negative_nodes = [int(i) for i, value in labels.items() if value < 0]
            if not negative_nodes:
                return current - rho_hat, current_w

            state = tuple(sorted((int(i), labels[i]) for i in labels.keys())) + ((-1, current.grade),)
            if state in visited:
                break
            visited.add(state)

            node = min(negative_nodes)
            current = current.simple_reflection(node)
            current_w = weyl_group.simple_reflection(node) * current_w

        raise ValueError("Failed to find dominant Lambda via affine simple reflections")

    @staticmethod
    def _apply_element_to_weight(
        algebra: "AffineLieAlgebra", element: Any, weight: "AffineWeight"
    ) -> "AffineWeight":
        return _apply_affine_element_to_weight(algebra, element, weight)

    @classmethod
    def _collect_quotient_representatives(
        cls,
        algebra: "AffineLieAlgebra",
        Lambda_hat: "AffineWeight",
        *,
        candidates: Iterable[Any],
    ) -> List[Any]:
        from .affine_weight import AffineWeight

        rho_hat = algebra.affine_rho()
        target = Lambda_hat + rho_hat
        target_domain = target.to_sagemath()
        by_weight: Dict[Tuple[Tuple[int, Any], ...], Any] = {}
        for w in candidates:
            try:
                acted = w.action(target_domain)
                image = AffineWeight.from_sagemath(algebra, acted) - rho_hat
            except Exception:
                image = cls._apply_element_to_weight(algebra, w, target) - rho_hat
            key = tuple(sorted(image.dynkin_labels().items())) + ((-1, image.grade),)
            current = by_weight.get(key)
            if current is None or int(w.length()) < int(current.length()):
                by_weight[key] = w
        return sorted(by_weight.values(), key=lambda w: (int(w.length()), tuple(_element_word_list(w))))

    @classmethod
    def _collect_stabilizer_and_quotient_representatives(
        cls,
        algebra: "AffineLieAlgebra",
        Lambda_hat: "AffineWeight",
        *,
        candidates: Iterable[Any],
    ) -> Tuple[List[Any], List[Any]]:
        """Collect stabilizer candidates and quotient representatives in one pass.

        This preserves the previous semantics while avoiding duplicated scans of
        the same bounded affine candidate set.
        """
        rho_hat = algebra.affine_rho()
        target = Lambda_hat + rho_hat
        target_domain = target.to_sagemath()

        by_weight: Dict[Any, Any] = {}
        stabilizer: List[Any] = []

        for w in candidates:
            try:
                acted = w.action(target_domain)
                if acted == target_domain:
                    stabilizer.append(w)

                # Fast key path: avoid expensive AffineWeight conversion on each
                # bounded affine candidate.
                key = tuple(acted.to_vector())
            except Exception:
                image = cls._apply_element_to_weight(algebra, w, target) - rho_hat
                if image == Lambda_hat:
                    stabilizer.append(w)
                key = tuple(sorted(image.dynkin_labels().items())) + ((-1, image.grade),)

            current = by_weight.get(key)
            if current is None or int(w.length()) < int(current.length()):
                by_weight[key] = w

        if stabilizer:
            identity = stabilizer[0].parent().one()
        elif by_weight:
            identity = next(iter(by_weight.values())).parent().one()
        else:
            identity = algebra.affine_weyl_group_sage().one()

        identity_word = tuple(_element_word_list(identity))
        if all(tuple(_element_word_list(w)) != identity_word for w in stabilizer):
            stabilizer.append(identity)

        stabilizer_sorted = sorted(
            stabilizer,
            key=lambda w: (int(w.length()), tuple(_element_word_list(w))),
        )
        quotient_sorted = sorted(
            by_weight.values(),
            key=lambda w: (int(w.length()), tuple(_element_word_list(w))),
        )
        return stabilizer_sorted, quotient_sorted


    def _build_W_affine_as_words_direct(
        self,
        normalized_translation_vectors: Iterable[Any],
    ) -> List[Any]:
        if self._finite_affine_elements_cache is None:
            finite_words = [
                tuple(int(i) for i in w.reduced_word())
                for w in list(self.algebra.affine_weyl_group()._finite_weyl_group)
            ]
            self._finite_affine_elements_cache = [
                self.algebra.affine_weyl_group_sage().from_reduced_word(list(word))
                for word in finite_words
            ]

        return _build_W_affine_as_words_direct_impl(
            self.algebra,
            normalized_translation_vectors,
            finite_affine_elements_cache=self._finite_affine_elements_cache,
        )

    def _build_candidates_and_collect_data(
        self,
        coroots: Iterable[Any],
        Lambda_hat: "AffineWeight",
    ) -> Tuple[List[Any], List[Any], List[Any]]:
        """Build candidates and collect stabilizer/quotient data in one pass."""
        affine_weyl_group = self.algebra.affine_weyl_group()
        affine_weyl_group_sage = self.algebra.affine_weyl_group_sage()

        rho_hat = self.algebra.affine_rho()
        dominant_weight = Lambda_hat + rho_hat
        dominant_weight_sage = dominant_weight.to_sagemath()
        finite_simple_reflections = [
            affine_weyl_group_sage.simple_reflection(i) for i in range(1, self.algebra.rank + 1)
        ]

        candidates: List[Any] = []
        by_weight: Dict[Any, Any] = {}
        stabilizer: List[Any] = []

        for beta in coroots:
            translation_word = tuple(int(i) for i in affine_weyl_group.translation_word_list(beta))
            translation_sage = affine_weyl_group_sage.from_reduced_word(list(translation_word))
            translated_dominant_weight = translation_sage.action(dominant_weight_sage)

            visited: set[Any] = set()
            queue = deque([(translated_dominant_weight, affine_weyl_group_sage.one())])

            while queue:
                acted, finite_affine = queue.popleft()
                state_key = tuple(acted.to_vector())
                if state_key in visited:
                    continue
                visited.add(state_key)

                candidate = finite_affine * translation_sage
                candidates.append(candidate)

                if acted == dominant_weight_sage:
                    stabilizer.append(candidate)

                current = by_weight.get(state_key)
                if current is None:
                    by_weight[state_key] = candidate
                elif int(candidate.length()) < int(current.length()):
                    by_weight[state_key] = candidate

                for simple_reflection in finite_simple_reflections:
                    queue.append((simple_reflection.action(acted), simple_reflection * finite_affine))

        if stabilizer:
            identity = stabilizer[0].parent().one()
        elif by_weight:
            identity = next(iter(by_weight.values())).parent().one()
        else:
            identity = self.algebra.affine_weyl_group_sage().one()

        identity_word = tuple(_element_word_list(identity))
        if all(tuple(_element_word_list(w)) != identity_word for w in stabilizer):
            stabilizer.append(identity)

        stabilizer_sorted = sorted(
            stabilizer,
            key=lambda w: (int(w.length()), tuple(_element_word_list(w))),
        )
        quotient_representatives = list(by_weight.values())
        return candidates, stabilizer_sorted, quotient_representatives


    def prepare_data(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translations: Optional[Iterable[Any]] = None,
        manual_translations: Optional[Iterable[Any]] = None,
    ) -> KazhdanLusztigData:
        rho_hat = self.algebra.affine_rho()
        Lambda_hat, w_to_Lambda = self._find_dominant_Lambda(lambda_hat)
        w_to_lambda = w_to_Lambda.inverse()

        # Λ + ρ (= w.λ) should be deominant
        dominant_weight = Lambda_hat + rho_hat

        if translations is not None and manual_translations is not None:
            raise ValueError("Pass either translations or manual_translations, not both")

        translation_source = manual_translations if manual_translations is not None else translations

        # ``translation_source`` accepts mixed explicit inputs (translation elements, coroot vectors,
        # or 0/None). Normalization happens in the affine Weyl group helper so the rest of the
        # pipeline always sees canonical coroot-space vectors.
        translation_inputs = list(translation_source) if translation_source is not None else None
        affine_weyl_group = self.algebra.affine_weyl_group()

        
        if translation_inputs is None:
            translation_order = order + max(
                0,
                int(QQ(dominant_weight.grade) - QQ(lambda_hat.grade)),
            )
            max_neg_shift = QQ(order) + QQ(dominant_weight.grade)
            
            # _translations_by_n_shift 返回 list of AffineWeylGroupSemidirectElement
            translation_inputs = self._translations_by_n_shift(
                dominant_weight,
                order=translation_order,
                max_neg_shift=max_neg_shift,
            )

        # 获取 translations 对应的 list of coroots
        coroots = affine_weyl_group._translations_to_coroots(
            translations=translation_inputs,
        )

        normalized_manual_translations = [
            affine_weyl_group.translation(beta) for beta in coroots
        ]

        (
            W_affine_as_words_sorted,
            stabilizer_candidates,
            quotient_representatives,
        ) = self._build_candidates_and_collect_data(
            coroots,
            Lambda_hat,
        )

        return KazhdanLusztigData(
            algebra=self.algebra,
            lambda_hat=lambda_hat,
            Lambda_hat=Lambda_hat,
            w_to_lambda=w_to_lambda,
            order=order,
            translations=normalized_manual_translations,
            W_affine_as_words=W_affine_as_words_sorted,
            stabilizer_candidates=stabilizer_candidates,
            quotient_representatives=quotient_representatives,
        )

    def numerator_terms(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translations: Optional[Iterable[Any]] = None,
        manual_translations: Optional[Iterable[Any]] = None,
    ) -> list[Any]:
        context = self.prepare_data(
            lambda_hat,
            order=order,
            translations=translations,
            manual_translations=manual_translations,
        )
        terms: list[KLNumeratorTerm] = []
        lower = context.w_to_lambda
        for representative in context.weyl_to_be_summed():
            coefficient = self.kl.affine_bounded_parabolic_Q_tilde(
                lower,
                representative,
                candidates=context.W_affine_as_words,
                stabilizer_candidates=context.stabilizer_candidates,
                at_one=True,
            )
            if coefficient == 0:
                continue
            weight = context.quotient_weight(representative)
            terms.append(
                KLNumeratorTerm(
                    representative=representative,
                    weight=weight,
                    coefficient=coefficient,
                )
            )
        return terms

    def character(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translations: Optional[Iterable[Any]] = None,
        manual_translations: Optional[Iterable[Any]] = None,
    ) -> FormalCharacter:
        result = FormalCharacter({}, max_grade=order, algebra=self.algebra)
        for term in self.numerator_terms(
            lambda_hat,
            order=order,
            translations=translations,
            manual_translations=manual_translations,
        ):
            verma = VermaCharacter(self.algebra, term.weight)
            result = result + term.coefficient * verma.character(max_grade=order)
        return result.truncate(order)

    def character_as_weight_space(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translations: Optional[Iterable[Any]] = None,
        manual_translations: Optional[Iterable[Any]] = None,
    ) -> Any:
        """Compute character as weight space decomposition (q-series with z variables).
        
        Returns the same format as IntegrableModuleCharacter.character():
        a symbolic expression in q, z1, z2, ... representing the weight space decomposition.
        
        This method converts the KL formula result (Verma module linear combination)
        into the weight space decomposition by using the IntegrableModuleCharacter
        machinery.
        """
        int_char = IntegrableModuleCharacter(lambda_hat)
        
        if manual_translations is not None:
            return int_char.character(order, manual_translations=manual_translations)
        
        return int_char.character(order)


def character_from_weight(
    algebra: "AffineLieAlgebra",
    weight: "AffineWeight",
    max_grade: int = 10,
) -> FormalCharacter:
    return VermaCharacter(algebra, weight).character(max_grade)
