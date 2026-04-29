from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import product
import time
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Tuple

from sage.all import IntegrableRepresentation as SageIntegrableRepresentation
from sage.all import Integer, QQ, SR, ZZ, binomial, matrix, prod, var, vector

if TYPE_CHECKING:
    from .affine_lie_algebra import AffineLieAlgebra
    from .affine_weight import AffineWeight


def _get_progress_bar(enabled: bool = False) -> Any:
    """Return a tqdm-compatible progress-bar wrapper or a no-op identity.

    When *enabled* is False (the default) returns a plain passthrough
    ``lambda it, **_: it`` so that no dependency on ``tqdm`` is
    triggered unless explicitly requested.
    """
    if not enabled:
        return lambda it, **_: it
    try:
        from tqdm import tqdm  # type: ignore[import-not-found]

        return tqdm
    except ImportError:
        return lambda it, **_: it


def _debug_log(enabled: bool, *args: Any) -> None:
    """Print a debug message when *enabled* is True."""
    if enabled:
        print(*args)


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
            rebuilt_highest_weight += (
                integer_coefficient * extended_affine_fundamental_weights[index]
            )
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


def _affine_dynkin_finite_labels(
    algebra: "AffineLieAlgebra", weight: "AffineWeight"
) -> List[Any]:
    sage_weight = weight.to_sagemath(extended=False)
    affine_dynkin_entries = list(sage_weight.to_vector()[0 : algebra.rank + 1])
    return [QQ(affine_dynkin_entries[i]) for i in range(1, algebra.rank + 1)]


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
        max_neg_shift_value = QQ(order)
    else:
        max_neg_shift_value = QQ(max_neg_shift)
    if max_neg_shift_value < 0:
        return {"translations": [], "stats": {}} if return_stats else []

    level = QQ(weight.level)
    if level == 0:
        raise ValueError("Translation enumeration by n-shift requires non-zero level")

    linear_coeffs = _affine_dynkin_finite_labels(algebra, weight)
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
        return {"translations": [t0], "stats": {"radius": 0}} if return_stats else [t0]

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

            current_const = (
                old_const
                + current_b[depth] * value_qq
                + level * gram_qq[depth][depth] * value_qq * value_qq / QQ(2)
            )
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


# NOTE: _translation_neg_shift 计算从 weight 出发，哪些 translations
# 让 n(weight) (也就是 weight.grade) 变化量 0 <= - Δn <= 上界
# _translation_neg_shift 只应该接收一个上界参数
# 现在 _translation_neg_shift 接受两个上界参数 order, max_neg_shift，显然有问题
# 在 KL 计算中，需要计算 0 <= - Δn <= order + n(Λhat + ρhat) - n(λhat)
# 因此 上界 参数应该由 order + n(Λhat + ρhat) - n(λhat) 给出
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
        max_neg_shift_value = QQ(order)
    else:
        max_neg_shift_value = QQ(max_neg_shift)
    if max_neg_shift_value < 0:
        return []

    if QQ(weight.level) == 0:
        raise ValueError("Translation enumeration by n-shift requires non-zero level")

    linear_coeffs = _affine_dynkin_finite_labels(algebra, weight)
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

    translations = [selected[key] for key in sorted(selected.keys())]
    print("translations = %s" % ([t.reduced_word() for t in translations]))
    return translations


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
            affine_weyl_group_sage.from_reduced_word(list(word)) for word in finite_words
        ]

    candidates: List[Any] = []
    for beta in normalized_translation_vectors:
        translation_word = tuple(int(i) for i in affine_weyl_group.translation_word_list(beta))
        translation_sage = affine_weyl_group_sage.from_reduced_word(list(translation_word))
        for finite_affine in finite_affine_elements:
            candidates.append(finite_affine * translation_sage)

    return candidates


def _sage_weight_to_affine(
    alg: "AffineLieAlgebra", sage_wt: Any, grade: Optional[Any] = None
) -> "AffineWeight":
    from .affine_weight import AffineWeight

    inferred_grade = QQ(grade) if grade is not None else QQ(0)
    try:
        vector_entries = list(sage_wt.to_vector())
        if grade is None and vector_entries:
            inferred_grade = QQ(vector_entries[-1])
    except Exception:
        pass

    return AffineWeight.from_sagemath(alg, sage_wt, grade=inferred_grade)


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
    w_to_Lambda: Optional[Any] = None

    @property
    def manual_translations(self) -> List[Any]:
        """Backward-compatible alias for older call sites/tests."""
        return self.translations

    @property
    def WLambda0(self) -> List[Any]:
        """Backward-compatible legacy alias for the affine stabilizer list."""
        return self.stabilizer_candidates

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
        if not self.quotient_representatives:
            return []
        lower = self.w_to_lambda
        result = [rep for rep in self.quotient_representatives if lower.bruhat_le(rep)]
        return sorted(result, key=lambda w: (int(w.length()), tuple(_element_word_list(w))))


class IntegrableModuleCharacter:
    """
    Compute character of integrable module of affine Lie algebra.

    Used ONLY when the affine highest weight has **non-negative** dynkin labels
    i.e.,  [λ0>=0, λ1>=0, ..., λr>=0]
    """

    def __init__(self, algebra: "AffineLieAlgebra") -> None:
        from .affine_lie_algebra import AffineLieAlgebra

        if not isinstance(algebra, AffineLieAlgebra):
            raise TypeError("IntegrableModuleCharacter expects an AffineLieAlgebra")
        self.algebra = algebra

    def _normalize_highest_weight(self, highest_weight: Any) -> tuple["AffineWeight", Any]:
        from .affine_weight import AffineWeight

        if isinstance(highest_weight, AffineWeight):
            if tuple(highest_weight.algebra._cartan_type) != tuple(self.algebra._cartan_type):
                raise ValueError(
                    "highest_weight algebra does not match IntegrableModuleCharacter algebra"
                )
            sage_weight = _to_extended_affine_highest_weight(
                self.algebra,
                highest_weight.to_sagemath(),
            )
            affine_weight = AffineWeight.from_sagemath(self.algebra, sage_weight, grade=0)
            return affine_weight, sage_weight

        cartan_type = list(highest_weight.parent().cartan_type())
        if tuple(cartan_type) != tuple(self.algebra._cartan_type):
            raise ValueError(
                "highest_weight algebra does not match IntegrableModuleCharacter algebra"
            )

        sage_weight = _to_extended_affine_highest_weight(self.algebra, highest_weight)
        affine_weight = AffineWeight.from_sagemath(self.algebra, sage_weight, grade=0)
        return affine_weight, sage_weight

    def _representation_for(self, highest_weight: Any) -> tuple["AffineWeight", Any]:
        highest_weight_affine, sage_weight = self._normalize_highest_weight(highest_weight)
        return highest_weight_affine, SageIntegrableRepresentation(sage_weight)

    def dominant_maximal_weights(self, highest_weight: Any) -> list[Any]:
        _, representation = self._representation_for(highest_weight)
        return list(representation.dominant_maximal_weights())

    def strings(self, highest_weight: Any, depth: int = 12) -> dict[Any, list[Any]]:
        d = int(depth)
        if d <= 0:
            raise ValueError("depth must be a positive integer")
        _, representation = self._representation_for(highest_weight)
        raw = representation.strings(d)
        return {weight: list(values) for weight, values in raw.items()}

    def from_weight(self, highest_weight: Any, weight: Any) -> tuple[int, ...]:
        _, representation = self._representation_for(highest_weight)
        return tuple(representation.from_weight(weight))

    def _sage_weight_to_affine_weight(self, highest_weight: Any, weight: Any) -> "AffineWeight":
        from .affine_weight import AffineWeight

        # In Sage's integrable highest-weight module coordinates, only alpha_0
        # contributes a delta-shift, so the grade is minus the alpha_0 count.
        root_coordinates = self.from_weight(highest_weight, weight)
        grade = -QQ(root_coordinates[0]) if root_coordinates else QQ(0)
        return AffineWeight.from_sagemath(self.algebra, weight, grade=grade)

    def _stable_strings_prefix(
        self,
        highest_weight: Any,
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
            data = self.strings(highest_weight, depth)
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
        return self.strings(highest_weight, depth)

    def _auto_translations(
        self, highest_weight: Any, dominant_weights: list[Any], *, order: int
    ) -> list[Any]:
        affine_weyl_group = self.algebra.affine_weyl_group()
        translation_vectors: dict[Tuple[int, ...], Any] = {}

        for weight in dominant_weights:
            affine_weight = self._sage_weight_to_affine_weight(highest_weight, weight)
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

        return [
            affine_weyl_group.translation(beta) for _, beta in sorted(translation_vectors.items())
        ]

    def _weyl_candidates(self, translations: Iterable[Any]) -> list[Any]:
        affine_weyl_group = self.algebra.affine_weyl_group()
        vectors = affine_weyl_group._translations_to_coroots(translations=translations)
        return _build_W_affine_as_words_direct_impl(self.algebra, vectors)

    def _orbit_representatives_for_weight(
        self,
        highest_weight: Any,
        weight: Any,
        *,
        candidates: Iterable[Any],
        order: int,
    ) -> list[tuple[Any, int]]:
        base_weight = self._sage_weight_to_affine_weight(highest_weight, weight)
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

        return sorted(
            selected, key=lambda item: (int(item[0].length()), tuple(_element_word_list(item[0])))
        )

    def _character_contribution(
        self,
        highest_weight: Any,
        weight: Any,
        representative: Any,
        string: list[Any],
        *,
        nmax: int,
        q_var: Any,
        z_vars: dict[int, Any],
    ) -> Any:
        from .affine_weight import AffineWeight

        base_weight = self._sage_weight_to_affine_weight(highest_weight, weight)
        delta = AffineWeight.delta(self.algebra)
        simple_roots = {
            i: AffineWeight.affine_simple_root(self.algebra, i)
            for i in range(1, self.algebra.rank + 1)
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

    def character(
        self,
        highest_weight: Any,
        order: int,
        *,
        show_progress: bool = False,
        debug: bool = False,
    ) -> Any:
        target_order = int(order)
        if target_order < 0:
            return 0

        tqdm_bar = _get_progress_bar(enabled=show_progress)
        logger = lambda *args: _debug_log(debug, *args)

        dominant_weights = self.dominant_maximal_weights(highest_weight)
        logger(
            f"[IntegrableModuleCharacter] Found {len(dominant_weights)} dominant maximal weight(s)"
        )

        logger("[IntegrableModuleCharacter] Computing auto translations …")
        translations = self._auto_translations(highest_weight, dominant_weights, order=target_order)
        logger(f"[IntegrableModuleCharacter] → {len(translations)} translation(s)")

        logger("[IntegrableModuleCharacter] Building Weyl-group candidates …")
        candidates = self._weyl_candidates(translations)
        logger(f"[IntegrableModuleCharacter] → {len(candidates)} candidate(s)")

        reps_by_weight: dict[Any, list[tuple[Any, int]]] = {}
        max_needed_depth_by_weight: dict[Any, int] = {}

        logger(
            "[IntegrableModuleCharacter] Computing orbit representatives for each dominant weight …"
        )
        for weight in dominant_weights:
            reps = self._orbit_representatives_for_weight(
                highest_weight,
                weight,
                candidates=candidates,
                order=target_order,
            )
            reps_by_weight[weight] = reps
            max_needed_depth_by_weight[weight] = max([nmax for _, nmax in reps], default=-1)

        if max(max_needed_depth_by_weight.values(), default=-1) < 0:
            return 0

        logger("[IntegrableModuleCharacter] Stabilizing string prefixes …")
        strings_data = self._stable_strings_prefix(highest_weight, max_needed_depth_by_weight)
        q_var = var("q")
        z_vars = {i: var(f"z{i}") for i in range(1, self.algebra.rank + 1)}

        total = SR(0)
        relevant_weights = [w for w in dominant_weights if max_needed_depth_by_weight[w] >= 0]
        logger(
            f"[IntegrableModuleCharacter] Summing contributions for {len(relevant_weights)} weight(s) …"
        )

        for weight in tqdm_bar(relevant_weights, desc="dominant weights", leave=False):
            string = strings_data[weight]
            for representative, nmax in tqdm_bar(
                reps_by_weight[weight],
                desc=f"reps (weight {weight})",
                leave=False,
            ):
                total += self._character_contribution(
                    highest_weight,
                    weight,
                    representative,
                    string,
                    nmax=nmax,
                    q_var=q_var,
                    z_vars=z_vars,
                )

        return sum(total.coefficient(q_var, n) * q_var**n for n in range(0, target_order + 1))

    def __repr__(self) -> str:
        return f"IntegrableModuleCharacter({self.algebra!r})"


class KazhdanLusztigCharacter:
    def __init__(self, algebra: "AffineLieAlgebra") -> None:
        from .kazhdan_lusztig import KazhdanLusztigPolynomials

        self.algebra = algebra
        self.kl_polynomial = KazhdanLusztigPolynomials(algebra.affine_weyl_group_sage())
        self._finite_affine_elements_cache: Optional[List[Any]] = None

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

    def _find_dominant_Lambda(
        self,
        lambda_hat: "AffineWeight",
        *,
        max_steps: int = 1000,
    ) -> Tuple["AffineWeight", Any]:
        rho_hat = self.algebra.affine_rho()
        weyl_group = self.kl_polynomial.weyl_group

        def legacy_affine_prefix_coefficients(weight: "AffineWeight") -> list[Any]:
            sage_weight = weight.to_sagemath(extended=False)
            return list(sage_weight.to_vector()[0 : self.algebra.rank])

        finite_coefficients = legacy_affine_prefix_coefficients(lambda_hat)
        if all(coeff > 0 for coeff in finite_coefficients):
            identity = weyl_group.one()
            return lambda_hat, identity

        lambda_plus_rho = lambda_hat + rho_hat
        checked = 0
        max_length = 0

        while checked < max_steps:
            elements_of_length = list(weyl_group.elements_of_length(max_length))
            elements_of_length.sort(key=lambda w: tuple(int(i) for i in w.reduced_word()))
            for w_to_Lambda in elements_of_length:
                checked += 1
                acted_weight = (
                    self._apply_element_to_weight(
                        self.algebra,
                        w_to_Lambda,
                        lambda_plus_rho,
                    )
                    - rho_hat
                )
                reduced_coefficients = legacy_affine_prefix_coefficients(acted_weight)
                if all(coeff >= -1 for coeff in reduced_coefficients):
                    return acted_weight, w_to_Lambda
                if checked >= max_steps:
                    break
            max_length += 1

        raise ValueError("Failed to find dominant Lambda by bounded affine Weyl search")

    @staticmethod
    def _apply_element_to_weight(
        algebra: "AffineLieAlgebra", element: Any, weight: "AffineWeight"
    ) -> "AffineWeight":
        return _apply_affine_element_to_weight(algebra, element, weight)

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

    def _denominator_candidates(self, order: int) -> List[Any]:
        rho_hat = self.algebra.affine_rho()
        affine_weyl_group = self.algebra.affine_weyl_group()
        denominator_translations = _translations_by_n_shift_impl(
            self.algebra,
            rho_hat,
            order=order,
            max_neg_shift=order,
        )
        denominator_coroots = affine_weyl_group._translations_to_coroots(
            translations=denominator_translations,
        )
        return self._build_W_affine_as_words_direct(denominator_coroots)

    def character_weight_list(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        Lambda_hat: Optional["AffineWeight"] = None,
        w_to_lambda: Optional[Any] = None,
        translations: Optional[Iterable[Any]] = None,
        show_progress: bool = False,
        debug: bool = False,
    ) -> list[Any]:
        """Port of CharacterNum from demos/MyAlgebra.py.

        Strictly mirrors the 8-step business logic of CharacterNum
        without using prepare_data::

            1. Determine Λ and w_T^{-1}λ
            2. Compute translation set T
            3. Construct W = W_fin × T
            4. Compute stabilizer W_{Λ,0}
            5. Compute dot-orbit w·(Λ+ρ)-ρ
            6. Coset deduplication (keep shortest rep)
            7. Bruhat filter: w_T^{-1}λ ≤ w'
            8. Compute Q̃ coefficients

        Returns a list of ``{weight: coefficient}`` dicts,
        matching CharacterNum's output format exactly.
        """
        from .affine_weight import AffineWeight

        tqdm_bar = _get_progress_bar(enabled=show_progress)
        logger = lambda *args: _debug_log(debug, *args)

        # ── Step 1: Determine Λ and w_T^{-1}λ ──
        # Mirrors CharacterNum lines 960-969
        if Lambda_hat is not None and w_to_lambda is not None:
            pass
        else:
            Lambda_hat, w_to_Lambda_hat = self._find_dominant_Lambda(lambda_hat)
            w_to_lambda = w_to_Lambda_hat.inverse()

        rho_hat = self.algebra.affine_rho()
        Lambda_plus_rho = Lambda_hat + rho_hat

        logger(
            f"[character_weight_list] λ̂ = {lambda_hat},  Λ̂ = {Lambda_hat},  order = {order}"
        )
        print(
            f"[character_weight_list] Computing KL numerator: "
            f"λ̂ = {lambda_hat},  Λ̂ = {Lambda_hat},  order = {order}",
            flush=True,
        )

        # ── Step 2: Compute translations ──
        # Mirrors CharacterNum lines 970-987
        # order_base = n(Λ+ρ) - n(λ)
        order_base = int(QQ(Lambda_plus_rho.grade) - QQ(lambda_hat.grade))
        translation_order = order_base + order

        if translations is None:
            print(
                f"[character_weight_list] computing translations "
                f"(order_base={order_base}, translation_order={translation_order}) …",
                flush=True,
            )
            translation_elements = self._translations_by_n_shift(
                Lambda_plus_rho,
                order=translation_order,
                max_neg_shift=QQ(translation_order),
            )
            print(
                f"[character_weight_list] → {len(translation_elements)} translation(s)",
                flush=True,
            )
        else:
            translation_elements = list(translations)
            print(
                f"[character_weight_list] using {len(translation_elements)} provided translation(s)",
                flush=True,
            )

        # ── Step 3: Construct W = W_fin × T ──
        # Mirrors CharacterNum lines 988-991
        # Convert translation elements to coroot vectors, then build W_affine
        affine_weyl_group = self.algebra.affine_weyl_group()
        coroots = affine_weyl_group._translations_to_coroots(
            translations=translation_elements,
        )

        print("[character_weight_list] building finite×translation candidates …", flush=True)
        W_affine_as_words = self._build_W_affine_as_words_direct(coroots)
        print(
            f"[character_weight_list] → {len(W_affine_as_words)} candidate(s)",
            flush=True,
        )

        # ── Step 4: Compute stabilizer W_{Λ,0} ──
        # Mirrors CharacterNum lines 993-994
        print("[character_weight_list] computing stabilizer WΛ₀ …", flush=True)
        stabilizer, quotient_representatives = self._collect_stabilizer_and_quotient_representatives(
            self.algebra,
            Lambda_hat,
            candidates=W_affine_as_words,
        )
        print(
            f"[character_weight_list] → {len(stabilizer)} stabilizer(s), "
            f"{len(quotient_representatives)} quotient rep(s)",
            flush=True,
        )

        # ── Step 5: Compute dot-orbit ──
        # Mirrors CharacterNum lines 998-1005
        # dot-action: w·(Λ+ρ)-ρ for all w in W
        rho_sage = rho_hat.to_sagemath(extended=True)
        target_sage = Lambda_plus_rho.to_sagemath(extended=True)

        print("[character_weight_list] computing dot-orbit …", flush=True)
        lambda_orbit_under_weyl_dot = []
        for w in tqdm_bar(W_affine_as_words, desc="dot-orbit", leave=False):
            acted = w.action(target_sage) - rho_sage
            lambda_orbit_under_weyl_dot.append(acted)
        print(
            f"[character_weight_list] → {len(lambda_orbit_under_weyl_dot)} orbit element(s)",
            flush=True,
        )

        # ── Step 6: Coset deduplication ──
        # Mirrors CharacterNum lines 1007-1032
        # For each unique image weight, keep the shortest Weyl representative
        print("[character_weight_list] building cosets …", flush=True)
        weights_to_be_summed: list[Any] = []
        cosets: list[Any] = []
        weight_index_by_key: Dict[Tuple[Any, ...], int] = {}

        for i, weight in enumerate(lambda_orbit_under_weyl_dot):
            key = tuple(weight.to_vector())
            w = W_affine_as_words[i]
            current_index = weight_index_by_key.get(key)
            if current_index is None:
                weight_index_by_key[key] = len(weights_to_be_summed)
                weights_to_be_summed.append(weight)
                cosets.append(w)
            elif int(w.length()) < int(cosets[current_index].length()):
                cosets[current_index] = w

        print(
            f"[character_weight_list] → {len(weights_to_be_summed)} unique weight(s)",
            flush=True,
        )

        # ── Step 7: Bruhat filter ──
        # Mirrors CharacterNum lines 1034-1039
        # Keep only coset representatives where w_T^{-1}λ ≤ w'
        print("[character_weight_list] filtering by Bruhat order …", flush=True)
        weyl_to_be_summed = [wp for wp in cosets if w_to_lambda.bruhat_le(wp)]
        print(
            f"[character_weight_list] → {len(weyl_to_be_summed)} Weyl element(s) to sum",
            flush=True,
        )

        # ── Step 8: Compute Q̃ coefficients ──
        # Mirrors CharacterNum lines 1047-1068
        print("[character_weight_list] computing Q̃ coefficients …", flush=True)
        result: list[Any] = []

        for representative in tqdm_bar(weyl_to_be_summed, desc="Q̃ computation", leave=True):
            coefficient = self.kl_polynomial.Q_tilde(
                w_to_lambda,
                representative,
                stabilizer_candidates=stabilizer,
            )
            logger(
                f"[character_weight_list]   w = {representative}  "
                f"(len={representative.length()})  Q̃ = {coefficient}"
            )
            if coefficient == 0:
                continue

            # Compute the image weight: w'·(Λ+ρ)-ρ
            acted_weight = representative.action(target_sage) - rho_sage
            weight = _sage_weight_to_affine(self.algebra, acted_weight)

            result.append({weight: coefficient})

        print(
            f"[character_weight_list] → {len(result)} non-zero term(s)",
            flush=True,
        )
        return result

    def numerator_q_series(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translations: Optional[Iterable[Any]] = None,
        manual_translations: Optional[Iterable[Any]] = None,
        show_progress: bool = False,
        debug: bool = False,
    ) -> Any:
        if translations is not None and manual_translations is not None:
            raise ValueError("Pass either translations or manual_translations, not both")
        translation_source = (
            manual_translations if manual_translations is not None else translations
        )

        legacy_terms = self.character_weight_list(
            lambda_hat,
            order=order,
            translations=translation_source,
            show_progress=show_progress,
            debug=debug,
        )
        q = var("q")
        # DON'T FORGET TO USE .subs({q:1}) for the coefficients
        numerator = sum(
            [
                SR(coefficient).subs(
                    {q: 1}) * self._character_contribution_from_weight(weight)
                for entry in legacy_terms
                for weight, coefficient in entry.items()
            ]
        )
        print(
            f"[KL] numerator_q_series: built from {len(legacy_terms)} weight term(s)",
            flush=True,
        )
        return numerator

    def _character_contribution_from_weight(self, weight: "AffineWeight") -> Any:
        algebra = self.algebra
        finite_rank = algebra.rank
        weight_lattice = algebra.affine_weight_lattice_sage()
        simple_roots = weight_lattice.simple_roots()
        variables = [var(f"b{i}") for i in range(0, finite_rank + 1)]
        q = var("q")
        sage_weight = weight.to_sagemath(extended=True)
        return prod(
            [
                variables[i]
                ** algebra.scalar_product(
                    weight,
                    _sage_weight_to_affine(algebra, simple_roots[i]),
                )
                for i in range(1, finite_rank + 1)
            ]
        ) * q ** (-QQ(weight.grade))

    def denominator_weight_list(self, order: int) -> List[Dict["AffineWeight", Any]]:
        rho_hat = self.algebra.affine_rho()
        denominator_candidates = self._denominator_candidates(order)
        rho_sage = rho_hat.to_sagemath(extended=True)
        denominator_weight_list = [
            {
                _sage_weight_to_affine(
                    self.algebra,
                    w.action(rho_sage) - rho_sage,
                ): (-1) ** (w.length() % 2)
            }
            for w in denominator_candidates
        ]
        print(
            f"[KL] denominator_q_series: built from {len(denominator_candidates)} Weyl element(s)",
            flush=True,
        )
        return denominator_weight_list

    def denominator_q_series(self, order: int) -> Any:
        denominator_weight_list = self.denominator_weight_list(order)
        denominator = sum(
            [
                SR(list(entry.values())[0])
                * self._character_contribution_from_weight(list(entry.keys())[0])
                for entry in denominator_weight_list
            ]
        )
        return denominator

    def character(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translations: Optional[Iterable[Any]] = None,
        manual_translations: Optional[Iterable[Any]] = None,
        show_progress: bool = False,
        debug: bool = False,
    ) -> Any:
        from sage.all import simplify as sage_simplify

        total_started = time.perf_counter()
        print(f"[KL] character: START  λ̂ = {lambda_hat},  order = {order}", flush=True)
        q = var("q")

        self.kl_polynomial.reset_profile_stats()
        self.kl_polynomial.set_profiling(True)

        numerator_started = time.perf_counter()
        numerator = self.numerator_q_series(
            lambda_hat,
            order=order,
            translations=translations,
            manual_translations=manual_translations,
            show_progress=show_progress,
            debug=debug,
        )
        print(numerator)
        numerator_seconds = time.perf_counter() - numerator_started

        denominator_started = time.perf_counter()
        denominator = self.denominator_q_series(order)
        print(denominator)
        denominator_seconds = time.perf_counter() - denominator_started

        ratio_started = time.perf_counter()
        print("[KL] character: computing ratio and Taylor expansion …", flush=True)
        

        result = sage_simplify((numerator / denominator).taylor(q, 0, order))
        ratio_seconds = time.perf_counter() - ratio_started

        total_seconds = time.perf_counter() - total_started
        stats = self.kl_polynomial.profile_stats()
        q_calls = int(stats.get("Q_calls", 0))
        q_time = float(stats.get("Q_total_seconds", 0.0))
        invpol_calls = int(stats.get("Q_invpol_calls", 0))
        invpol_time = float(stats.get("Q_invpol_seconds", 0.0))
        q_tilde_calls = int(stats.get("Q_tilde_calls", 0))
        q_tilde_time = float(stats.get("Q_tilde_total_seconds", 0.0))
        q_tilde_terms = int(stats.get("Q_tilde_stabilizer_terms", 0))
        q_cache_hits = int(stats.get("Q_cache_hits_at_one", 0)) + int(
            stats.get("Q_cache_hits_poly", 0)
        )
        print(
            (
                "[KL][timing] total={:.3f}s | numerator={:.3f}s | denominator={:.3f}s "
                "| ratio+taylor={:.3f}s"
            ).format(total_seconds, numerator_seconds, denominator_seconds, ratio_seconds),
            flush=True,
        )
        print(
            (
                "[KL][timing] Q_tilde: calls={} terms={} time={:.3f}s | "
                "Q: calls={} cache_hits={} time={:.3f}s | "
                "invpol: calls={} time={:.3f}s"
            ).format(
                q_tilde_calls,
                q_tilde_terms,
                q_tilde_time,
                q_calls,
                q_cache_hits,
                q_time,
                invpol_calls,
                invpol_time,
            ),
            flush=True,
        )
        self.kl_polynomial.set_profiling(False)
        print("[KL] character: DONE", flush=True)
        return result
