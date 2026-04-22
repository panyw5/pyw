from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Tuple

from sage.all import Integer, QQ, RR, binomial, ceil, floor, matrix, vector
from sage.all import ZZ

if TYPE_CHECKING:
    from .affine_lie_algebra import AffineLieAlgebra
    from .affine_weight import AffineWeight


def _element_word_list(element: Any) -> list[int]:
    if hasattr(element, "reduced_word_list"):
        return list(element.reduced_word_list())
    return [int(i) for i in element.reduced_word()]


@dataclass(frozen=True)
class KLNumeratorTerm:
    representative: Any
    weight: "AffineWeight"
    coefficient: Any


@dataclass
class KazhdanLusztigData:
    algebra: Any
    lambda_hat: "AffineWeight"
    rho_hat: "AffineWeight"
    Lambda_hat: "AffineWeight"
    w_to_Lambda: Any
    w_to_lambda: Any
    order: int
    translation_bounds: Optional[Dict[int, Tuple[int, int]]]
    translations: List[Any]
    ambient_candidates: List[Any]
    integral_case: bool
    weyl_lambda_candidates: List[Any]
    stabilizer_candidates: List[Any]
    quotient_representatives: List[Any]

    def apply(self, element: Any, weight: "AffineWeight") -> "AffineWeight":
        semidirect = self.algebra.affine_weyl_group()
        word = tuple(_element_word_list(element))
        return semidirect.from_word(word).action(weight)

    def quotient_weight(self, representative: Any) -> "AffineWeight":
        return self.apply(representative, self.Lambda_hat + self.rho_hat) - self.rho_hat

    def weyl_to_be_summed(self) -> List[Any]:
        from .bruhat import BruhatOrder

        if not self.ambient_candidates:
            return []
        bruhat = BruhatOrder(self.ambient_candidates[0].parent())
        lower = self.w_to_lambda
        result = []
        for w in self.quotient_representatives:
            if bruhat.le(lower, w):
                result.append(w)
        return sorted(result, key=lambda w: (int(w.length()), tuple(_element_word_list(w))))


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
        self._algebra = algebra
        self._inverse_cache: Dict[int, FormalCharacter] = {}

    @property
    def algebra(self) -> "AffineLieAlgebra":
        return self._algebra

    def inverse(self, max_grade: int = 10) -> FormalCharacter:
        cached = self._inverse_cache.get(max_grade)
        if cached is not None:
            return cached

        inverse = self._compute_inverse_product(max_grade)
        self._inverse_cache[max_grade] = inverse
        return inverse

    def _compute_inverse_product(self, max_grade: int) -> FormalCharacter:
        rank = int(getattr(self._algebra, "rank", 1) or 1)
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

        return FormalCharacter(coeffs, max_grade=max_grade, algebra=self._algebra)


class VermaCharacter:
    def __init__(self, algebra: "AffineLieAlgebra", weight: "AffineWeight") -> None:
        self._algebra = algebra
        self._weight = weight
        self._denominator = WeylKacDenominator(algebra)

    @property
    def algebra(self) -> "AffineLieAlgebra":
        return self._algebra

    @property
    def weight(self) -> "AffineWeight":
        return self._weight

    def character(self, max_grade: int = 10) -> FormalCharacter:
        weight_grade = int(getattr(self._weight, "grade", 0))
        inverse = self._denominator.inverse(max_grade + abs(weight_grade))
        return inverse.shift(-weight_grade).truncate(max_grade)


class KazhdanLusztigCharacter:
    def __init__(self, algebra: "AffineLieAlgebra") -> None:
        from .kazhdan_lusztig import KazhdanLusztigPolynomials

        self._algebra = algebra
        self._kl_polynomial = KazhdanLusztigPolynomials(algebra.affine_weyl_group_sage())

    @property
    def algebra(self) -> "AffineLieAlgebra":
        return self._algebra

    @property
    def kl(self):
        return self._kl_polynomial

    @staticmethod
    def _default_translation_bounds(weight: "AffineWeight", *, order: int) -> Dict[int, Tuple[int, int]]:
        labels = weight.dynkin_labels()
        finite_indices = [int(i) for i in labels.keys() if int(i) != 0]
        bounds: Dict[int, Tuple[int, int]] = {}
        for i in finite_indices:
            label = labels.get(i, 0)
            lower = -int(order)
            upper = int(order + max(0, label))
            bounds[i] = (lower, upper)
        return bounds

    @staticmethod
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

    def _exact_translations_by_n_shift(
        self,
        weight: "AffineWeight",
        *,
        order: int,
    ) -> List[Any]:
        semidirect = self._algebra.affine_weyl_group()
        idxs = [int(i) for i in semidirect._finite_coroot_space.index_set()]

        if QQ(weight.level) <= 0:
            raise ValueError(
                "Exact translation selection requires positive level for Lambda_hat + rho_hat; "
                "pass translation_bounds explicitly to use a manual bounded search"
            )

        basis = semidirect._finite_coroot_space.simple_roots()
        labels = weight.dynkin_labels()
        linear = vector(QQ, [QQ(labels.get(i, 0)) for i in idxs])
        gram = matrix(
            QQ,
            len(idxs),
            len(idxs),
            lambda r, c: QQ(basis[idxs[r]].to_ambient().inner_product(basis[idxs[c]].to_ambient())),
        )
        quadratic = (QQ(weight.level) / 2) * gram

        if not quadratic.is_positive_definite():
            raise ValueError(
                "Exact translation selection requires a positive-definite quadratic form; "
                "pass translation_bounds explicitly to use a manual bounded search"
            )

        inverse = quadratic.inverse()
        center = inverse * linear / 2
        radius_sq = QQ(order) + (center.dot_product(quadratic * center))

        if radius_sq < 0:
            return []

        ranges: list[range] = []
        for i, center_i in enumerate(center):
            projection_sq = radius_sq * inverse[i, i]
            bound = QQ(self._ceil_sqrt_qq(projection_sq))
            lower = int(floor(-center_i - bound))
            upper = int(ceil(-center_i + bound))
            ranges.append(range(lower, upper + 1))

        max_neg_shift = QQ(order) + QQ(weight.grade)
        if max_neg_shift < 0:
            return []

        selected: dict[tuple[int, ...], Any] = {}
        for coeffs in product(*ranges):
            beta = semidirect._zero_beta
            for i, coeff in zip(idxs, coeffs):
                if coeff:
                    beta += int(coeff) * basis[i]

            neg_shift = self._translation_neg_shift(weight, beta)
            if neg_shift < 0 or neg_shift > max_neg_shift:
                continue
            selected[tuple(int(coeff) for coeff in coeffs)] = semidirect.translation(beta)

        return [selected[key] for key in sorted(selected.keys())]

    def _translation_neg_shift(self, weight: "AffineWeight", beta: Any) -> Any:
        translated = self._algebra.affine_weyl_group().translation(beta).action(weight)
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
    def _is_integral_affine_weight(weight: "AffineWeight") -> bool:
        labels = weight.dynkin_labels()
        return all(value in ZZ for value in labels.values())

    @staticmethod
    def _apply_element_to_weight(
        algebra: "AffineLieAlgebra", element: Any, weight: "AffineWeight"
    ) -> "AffineWeight":
        semidirect = algebra.affine_weyl_group()
        word = tuple(_element_word_list(element))
        return semidirect.from_word(word).action(weight)

    @classmethod
    def _collect_quotient_representatives(
        cls,
        algebra: "AffineLieAlgebra",
        Lambda_hat: "AffineWeight",
        *,
        rho_hat: "AffineWeight",
        candidates: Iterable[Any],
    ) -> List[Any]:
        by_weight: Dict[Tuple[Tuple[int, Any], ...], Any] = {}
        for w in candidates:
            image = cls._apply_element_to_weight(algebra, w, Lambda_hat + rho_hat) - rho_hat
            key = tuple(sorted(image.dynkin_labels().items())) + ((-1, image.grade),)
            current = by_weight.get(key)
            if current is None or int(w.length()) < int(current.length()):
                by_weight[key] = w
        return sorted(by_weight.values(), key=lambda w: (int(w.length()), tuple(_element_word_list(w))))

    def build_context(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translation_bounds: Optional[Dict[int, Tuple[int, int]]] = None,
        translations: Optional[Iterable[Any]] = None,
    ) -> KazhdanLusztigData:
        rho_hat = self._algebra.affine_rho()
        Lambda_hat, w_to_Lambda = self._find_dominant_Lambda(lambda_hat)
        w_to_lambda = w_to_Lambda.inverse()

        if translation_bounds is not None and translations is not None:
            raise ValueError("Provide either translation_bounds or translations, not both")

        selected_weight = Lambda_hat + rho_hat
        selected_bounds = dict(translation_bounds) if translation_bounds is not None else None
        raw_translations = list(translations) if translations is not None else None
        semidirect = self._algebra.affine_weyl_group()

        if raw_translations is None and selected_bounds is None:
            raw_translations = self._exact_translations_by_n_shift(selected_weight, order=order)

        normalized_translation_vectors = semidirect._translation_vectors_from_inputs(
            translation_bounds=selected_bounds,
            translations=raw_translations,
        )
        normalized_translations = [semidirect.translation(beta) for beta in normalized_translation_vectors]

        ambient_candidates = self._kl_polynomial.affine_bounded_elements(
            self._algebra,
            translations=normalized_translations,
            factor_order="st",
        )
        integral_case = self._is_integral_affine_weight(Lambda_hat)
        weyl_lambda_candidates = list(ambient_candidates)
        stabilizer_candidates = self._kl_polynomial.affine_stabilizer(
            Lambda_hat.to_sagemath(),
            rho_hat=rho_hat.to_sagemath(),
            candidates=ambient_candidates,
            algebra=self._algebra,
        )
        quotient_representatives = self._collect_quotient_representatives(
            self._algebra,
            Lambda_hat,
            rho_hat=rho_hat,
            candidates=weyl_lambda_candidates,
        )

        return KazhdanLusztigData(
            algebra=self._algebra,
            lambda_hat=lambda_hat,
            rho_hat=rho_hat,
            Lambda_hat=Lambda_hat,
            w_to_Lambda=w_to_Lambda,
            w_to_lambda=w_to_lambda,
            order=order,
            translation_bounds=selected_bounds,
            translations=list(normalized_translations),
            ambient_candidates=ambient_candidates,
            integral_case=integral_case,
            weyl_lambda_candidates=weyl_lambda_candidates,
            stabilizer_candidates=stabilizer_candidates,
            quotient_representatives=quotient_representatives,
        )

    def numerator_terms(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
        translation_bounds: Optional[Dict[int, Tuple[int, int]]] = None,
        translations: Optional[Iterable[Any]] = None,
    ) -> list[Any]:
        context = self.build_context(
            lambda_hat,
            order=order,
            translation_bounds=translation_bounds,
            translations=translations,
        )
        terms: list[KLNumeratorTerm] = []
        lower = context.w_to_lambda
        for representative in context.weyl_to_be_summed():
            coefficient = self._kl_polynomial.affine_bounded_parabolic_Q_tilde(
                lower,
                representative,
                candidates=context.ambient_candidates,
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
        translation_bounds: Optional[Dict[int, Tuple[int, int]]] = None,
        translations: Optional[Iterable[Any]] = None,
    ) -> FormalCharacter:
        result = FormalCharacter({}, max_grade=order, algebra=self._algebra)
        for term in self.numerator_terms(
            lambda_hat,
            order=order,
            translation_bounds=translation_bounds,
            translations=translations,
        ):
            verma = VermaCharacter(self._algebra, term.weight)
            result = result + term.coefficient * verma.character(max_grade=order)
        return result.truncate(order)


def character_from_weight(
    algebra: "AffineLieAlgebra",
    weight: "AffineWeight",
    max_grade: int = 10,
) -> FormalCharacter:
    return VermaCharacter(algebra, weight).character(max_grade)
