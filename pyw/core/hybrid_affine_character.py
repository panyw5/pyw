from __future__ import annotations

import logging
import time
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Callable, Mapping

from sage.all import QQ, ZZ, WeylCharacterRing

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .affine_lie_algebra import AffineLieAlgebra
    from .affine_weight import AffineWeight, AffineWeightKey
    from .character import BoundedKLOrbit


def _require_untwisted_affine(algebra: "AffineLieAlgebra") -> None:
    cartan_type = algebra._cartan_type_obj
    if not algebra.is_affine:
        raise ValueError("The hybrid affine character helpers require an affine type")
    if cartan_type is None or algebra._finite_root_system is None:
        raise ValueError("The affine root-system data is not initialized")
    if hasattr(cartan_type, "is_untwisted_affine") and not cartan_type.is_untwisted_affine():
        raise NotImplementedError("Affine simple-root coordinates require an untwisted affine type")


@dataclass(frozen=True)
class KLCandidate:
    """A bounded KL orbit term whose weight is above a target weight."""

    weyl_key: tuple[int, ...]
    representative: Any
    orbit_weight: "AffineWeight"
    difference_coordinates: tuple[int, ...]


@dataclass(frozen=True)
class BoundedAffinePositiveRoot:
    """An affine positive root encoded in the affine simple-root basis."""

    coordinates: tuple[int, ...]
    multiplicity: int
    is_imaginary: bool


@dataclass(frozen=True)
class BoundedWeightTree:
    """Exact affine weights grouped by nonnegative depth below a highest weight."""

    highest_weight: "AffineWeight"
    order: int
    by_depth: Mapping[int, tuple["AffineWeight", ...]]
    weight_by_key: Mapping["AffineWeightKey", "AffineWeight"]

    def weights_at_depth(self, depth: int) -> tuple["AffineWeight", ...]:
        return self.by_depth.get(depth, ())

    def contains_key(self, key: "AffineWeightKey") -> bool:
        return key in self.weight_by_key


@dataclass
class HybridCharacterStats:
    """Measurements for one hybrid affine-character calculation."""

    freudenthal_targets: int = 0
    degenerate_targets: int = 0
    candidate_query_hits: int = 0
    partition_nonzero_hits: int = 0
    q_tilde_requests: int = 0
    q_tilde_unique_evaluations: int = 0
    q_tilde_cache_hits: int = 0
    orbit_seconds: float = 0.0
    weight_tree_seconds: float = 0.0
    partition_seconds: float = 0.0
    multiplicity_seconds: float = 0.0
    total_seconds: float = 0.0
    kl_q_calls: int = 0
    kl_q_cache_hits: int = 0
    kl_q_tilde_calls: int = 0

    def as_dict(self) -> dict[str, Any]:
        return dict(vars(self))


def affine_simple_root_coordinates(
    algebra: "AffineLieAlgebra",
    difference: "AffineWeight",
) -> tuple[Any, ...]:
    """Express a level-zero affine weight difference in affine simple roots."""
    _require_untwisted_affine(algebra)
    if difference.algebra != algebra:
        raise ValueError("The affine weight difference belongs to a different algebra")
    if difference.level != 0:
        raise ValueError("An affine root-lattice difference must have level zero")

    coefficient_alpha_0 = QQ(difference.grade)
    finite_root_system = algebra._finite_root_system
    if finite_root_system is None:
        raise ValueError("The affine root-system data is not initialized")
    finite_root_coordinates = (
        difference.finite_part_in_simple_root_basis().monomial_coefficients()
    )
    return (coefficient_alpha_0,) + tuple(
        QQ(finite_root_coordinates.get(index, 0))
        + QQ(algebra.marks.get(index, 0)) * coefficient_alpha_0
        for index in finite_root_system.root_space().index_set()
    )


def is_in_affine_positive_root_cone(
    algebra: "AffineLieAlgebra",
    difference: "AffineWeight",
) -> bool:
    """Return whether a difference is a nonnegative integral sum of simple roots."""
    _require_untwisted_affine(algebra)
    if difference.algebra != algebra or difference.level != 0 or difference.grade < 0:
        return False
    coordinates = affine_simple_root_coordinates(algebra, difference)
    return all(coordinate in ZZ and coordinate >= 0 for coordinate in coordinates)


def bounded_kl_candidates_above(
    orbit: "BoundedKLOrbit",
    target_weight: "AffineWeight",
) -> tuple[KLCandidate, ...]:
    """Find bounded KL orbit weights above a target in the affine root cone."""
    from .affine_weight import affine_weight_key

    algebra = orbit.lambda_hat.algebra
    _require_untwisted_affine(algebra)
    if target_weight.algebra != orbit.lambda_hat.algebra:
        raise ValueError("The target weight belongs to a different algebra")

    lambda_minus_target_coordinates = affine_simple_root_coordinates(
        algebra,
        orbit.lambda_hat - target_weight,
    )
    candidates: list[KLCandidate] = []
    for weyl_key, orbit_weight in orbit.weyl_to_weight.items():
        representative = orbit.weyl_to_representative.get(weyl_key)
        if representative is None:
            representative = orbit.weight_to_weyl[affine_weight_key(orbit_weight)]
        if orbit_weight.grade < target_weight.grade:
            continue

        relative_coordinates = orbit.weyl_to_relative_root_coordinates.get(weyl_key)
        if relative_coordinates is None:
            relative_coordinates = affine_simple_root_coordinates(
                algebra,
                orbit_weight - orbit.lambda_hat,
            )
        coordinates = tuple(
            relative_coordinate + target_coordinate
            for relative_coordinate, target_coordinate in zip(
                relative_coordinates,
                lambda_minus_target_coordinates,
                strict=True,
            )
        )
        if any(coordinate not in ZZ or coordinate < 0 for coordinate in coordinates):
            continue
        candidates.append(
            KLCandidate(
                weyl_key=weyl_key,
                representative=representative,
                orbit_weight=orbit_weight,
                difference_coordinates=tuple(int(coordinate) for coordinate in coordinates),
            )
        )
    return tuple(candidates)


def bounded_affine_positive_roots(
    algebra: "AffineLieAlgebra",
    max_grade: int,
) -> tuple[BoundedAffinePositiveRoot, ...]:
    """Enumerate untwisted affine positive roots through ``max_grade``."""
    _require_untwisted_affine(algebra)
    if max_grade not in ZZ or max_grade < 0:
        raise ValueError("The maximum affine grade must be a nonnegative integer")
    max_grade = int(max_grade)

    finite_root_system = algebra._finite_root_system
    if finite_root_system is None:
        raise ValueError("The affine root-system data is not initialized")
    finite_root_lattice = finite_root_system.root_lattice()
    finite_indices = tuple(finite_root_lattice.index_set())
    positive_finite_roots = tuple(finite_root_lattice.positive_roots())
    all_finite_roots = positive_finite_roots + tuple(-root for root in positive_finite_roots)
    marks = tuple(int(algebra.marks[index]) for index in finite_indices)

    def coordinates(finite_root: Any, grade: int) -> tuple[int, ...]:
        coefficients = finite_root.monomial_coefficients()
        return (grade,) + tuple(
            int(coefficients.get(index, 0)) + grade * mark
            for index, mark in zip(finite_indices, marks, strict=True)
        )

    roots = [
        BoundedAffinePositiveRoot(
            coordinates=coordinates(root, 0),
            multiplicity=1,
            is_imaginary=False,
        )
        for root in positive_finite_roots
    ]
    for grade in range(1, max_grade + 1):
        roots.extend(
            BoundedAffinePositiveRoot(
                coordinates=coordinates(root, grade),
                multiplicity=1,
                is_imaginary=False,
            )
            for root in all_finite_roots
        )
        roots.append(
            BoundedAffinePositiveRoot(
                coordinates=(grade,) + tuple(grade * mark for mark in marks),
                multiplicity=algebra.rank,
                is_imaginary=True,
            )
        )
    return tuple(roots)


class BoundedAffineVermaPartition:
    """Affine Kostant partition multiplicities inside a coordinate box."""

    def __init__(
        self,
        algebra: "AffineLieAlgebra",
        max_coordinates: Iterable[int],
    ) -> None:
        _require_untwisted_affine(algebra)
        raw_coordinates = tuple(max_coordinates)
        if any(coordinate not in ZZ for coordinate in raw_coordinates):
            raise ValueError("Affine simple-root coordinate bounds must be integers")
        coordinates = tuple(int(coordinate) for coordinate in raw_coordinates)
        if len(coordinates) != algebra.rank + 1:
            raise ValueError("The coordinate bound must have one entry per affine simple root")
        if any(coordinate < 0 for coordinate in coordinates):
            raise ValueError("Affine simple-root coordinate bounds must be nonnegative")

        self.algebra = algebra
        self.max_coordinates = coordinates
        self.positive_roots = tuple(
            root
            for root in bounded_affine_positive_roots(algebra, coordinates[0])
            if all(
                root_coordinate <= bound
                for root_coordinate, bound in zip(root.coordinates, coordinates, strict=True)
            )
        )
        self._multiplicities = self._build_table()

    def _build_table(self) -> dict[tuple[int, ...], int]:
        zero = (0,) * len(self.max_coordinates)
        multiplicities = {zero: 1}
        for root in self.positive_roots:
            for _ in range(root.multiplicity):
                updated = dict(multiplicities)
                for coordinates, coefficient in multiplicities.items():
                    shifted = tuple(
                        coordinate + root_coordinate
                        for coordinate, root_coordinate in zip(
                            coordinates, root.coordinates, strict=True
                        )
                    )
                    while all(
                        coordinate <= bound
                        for coordinate, bound in zip(
                            shifted, self.max_coordinates, strict=True
                        )
                    ):
                        updated[shifted] = updated.get(shifted, 0) + coefficient
                        shifted = tuple(
                            coordinate + root_coordinate
                            for coordinate, root_coordinate in zip(
                                shifted, root.coordinates, strict=True
                            )
                        )
                multiplicities = updated
            logger.debug(
                "Applied affine root coordinates=%s multiplicity=%s; table_size=%s",
                root.coordinates,
                root.multiplicity,
                len(multiplicities),
            )
        return multiplicities

    def multiplicity(self, coordinates: Iterable[int]) -> int:
        """Return the partition multiplicity, or zero outside the bounded cone."""
        key = tuple(coordinates)
        if len(key) != len(self.max_coordinates):
            raise ValueError("The query must have one coordinate per affine simple root")
        if any(coordinate not in ZZ or coordinate < 0 for coordinate in key):
            return 0
        integer_key = tuple(int(coordinate) for coordinate in key)
        if any(
            coordinate > bound
            for coordinate, bound in zip(
                integer_key, self.max_coordinates, strict=True
            )
        ):
            return 0
        return self._multiplicities.get(integer_key, 0)


def kl_coefficient_at_weight(
    target_weight: "AffineWeight",
    orbit: "BoundedKLOrbit",
    partition: BoundedAffineVermaPartition,
    kl_polynomial: Any,
    *,
    q_tilde_cache: dict[tuple[int, ...], Any],
    stats: HybridCharacterStats | None = None,
    candidates: tuple[KLCandidate, ...] | None = None,
) -> Any:
    """Extract one degenerate multiplicity with filtered, lazy Q-tilde calls."""
    if target_weight.algebra != orbit.lambda_hat.algebra:
        raise ValueError("The target weight belongs to a different algebra")
    if partition.algebra != orbit.lambda_hat.algebra:
        raise ValueError("The Verma partition belongs to a different algebra")

    value = QQ(0)
    if candidates is None:
        candidates = orbit.candidates_above(target_weight)
    if stats is not None:
        stats.candidate_query_hits += len(candidates)
    for candidate in candidates:
        partition_multiplicity = partition.multiplicity(candidate.difference_coordinates)
        if partition_multiplicity == 0:
            continue
        if stats is not None:
            stats.partition_nonzero_hits += 1
            stats.q_tilde_requests += 1

        if candidate.weyl_key in q_tilde_cache:
            coefficient = q_tilde_cache[candidate.weyl_key]
            if stats is not None:
                stats.q_tilde_cache_hits += 1
        else:
            coefficient = kl_polynomial.Q_tilde(
                orbit.w_to_lambda,
                candidate.representative,
                stabilizer_candidates=orbit.stabilizer,
            )
            q_tilde_cache[candidate.weyl_key] = coefficient
            if stats is not None:
                stats.q_tilde_unique_evaluations += 1
        value += QQ(coefficient) * partition_multiplicity
    return value


def build_bounded_weight_tree(
    highest_weight: "AffineWeight",
    *,
    order: int,
) -> BoundedWeightTree:
    """Build a bounded affine weight domain for a finite-integrable highest weight."""
    from .affine_weight import affine_weight_key

    algebra = highest_weight.algebra
    _require_untwisted_affine(algebra)
    if order not in ZZ or order < 0:
        raise ValueError("The character order must be a nonnegative integer")
    finite_root_system = algebra._finite_root_system
    if finite_root_system is None:
        raise ValueError("The affine root-system data is not initialized")
    highest_labels = highest_weight.dynkin_labels()
    finite_indices = tuple(finite_root_system.index_set())
    if any(
        highest_labels[index] not in ZZ or highest_labels[index] < 0
        for index in finite_indices
    ):
        raise NotImplementedError(
            "The bounded weight tree requires a highest weight that is "
            "dominant integral for the finite simple roots"
        )

    simple_roots = algebra.affine_simple_roots()
    finite_indices = tuple(finite_root_system.index_set())
    by_depth: dict[int, tuple["AffineWeight", ...]] = {}
    all_weights: dict["AffineWeightKey", "AffineWeight"] = {}
    previous_grade = (highest_weight,)

    for depth in range(int(order) + 1):
        seeds = (
            previous_grade
            if depth == 0
            else tuple(weight - simple_roots[0] for weight in previous_grade)
        )
        weights_at_depth = {affine_weight_key(weight): weight for weight in seeds}
        pending = list(seeds)
        while pending:
            weight = pending.pop()
            labels = weight.dynkin_labels()
            for index in finite_indices:
                label = labels[index]
                if label not in ZZ:
                    raise ValueError("Finite Dynkin labels must be integral")
                for step in range(1, max(0, int(label)) + 1):
                    descendant = weight - step * simple_roots[index]
                    key = affine_weight_key(descendant)
                    if key not in weights_at_depth:
                        weights_at_depth[key] = descendant
                        pending.append(descendant)

        ordered_weights = tuple(
            weights_at_depth[key]
            for key in sorted(weights_at_depth, key=repr)
        )
        by_depth[depth] = ordered_weights
        all_weights.update(weights_at_depth)
        previous_grade = ordered_weights

    return BoundedWeightTree(
        highest_weight=highest_weight,
        order=int(order),
        by_depth=MappingProxyType(by_depth),
        weight_by_key=MappingProxyType(all_weights),
    )


def freudenthal_denominator(
    highest_weight: "AffineWeight",
    target_weight: "AffineWeight",
) -> Any:
    """Return the exact shifted-norm difference in Freudenthal recursion."""
    from .affine_weight import AffineWeight

    if target_weight.algebra != highest_weight.algebra:
        raise ValueError("The target weight belongs to a different algebra")
    if target_weight.level != highest_weight.level:
        raise ValueError("The target and highest weights must have the same level")
    rho_hat = AffineWeight.rho_hat(highest_weight.algebra)
    return (highest_weight + rho_hat).norm_squared() - (
        target_weight + rho_hat
    ).norm_squared()


def _affine_root_from_coordinates(
    algebra: "AffineLieAlgebra",
    coordinates: tuple[int, ...],
) -> "AffineWeight":
    from .affine_weight import AffineWeight

    root = AffineWeight.zero(algebra)
    simple_roots = algebra.affine_simple_roots()
    for index, coefficient in enumerate(coordinates):
        root += coefficient * simple_roots[index]
    return root


def classically_dominant_affine_weight(target_weight: "AffineWeight") -> "AffineWeight":
    """Return the finite-Weyl dominant representative at fixed level and grade."""
    from .affine_weight import AffineWeight

    _require_untwisted_affine(target_weight.algebra)
    finite_part = target_weight.finite_part.to_dominant_chamber()
    return AffineWeight(
        target_weight.algebra,
        finite_part,
        target_weight.level,
        target_weight.grade,
    )


def decompose_finite_irreducibles(
    weight_multiplicities: Mapping["AffineWeight", int],
) -> dict["AffineWeight", int]:
    """Decompose one finite-Weyl-invariant affine weight space into irreducibles."""
    from .affine_weight import AffineWeight, affine_weight_key

    if not weight_multiplicities:
        return {}

    weights = tuple(weight_multiplicities)
    reference_weight = weights[0]
    algebra = reference_weight.algebra
    _require_untwisted_affine(algebra)
    if any(weight.algebra != algebra for weight in weights):
        raise ValueError("All weights must belong to the same affine algebra")
    if any(
        weight.level != reference_weight.level or weight.grade != reference_weight.grade
        for weight in weights
    ):
        raise ValueError("A finite decomposition requires one fixed affine level and grade")

    remaining = {
        affine_weight_key(weight): int(multiplicity)
        for weight, multiplicity in weight_multiplicities.items()
        if multiplicity != 0
    }
    if any(multiplicity < 0 for multiplicity in remaining.values()):
        raise ValueError("Weight multiplicities must be nonnegative")
    weights_by_key = {affine_weight_key(weight): weight for weight in weights}
    finite_type = algebra._finite_type
    finite_root_system = algebra._finite_root_system
    if finite_type is None or finite_root_system is None:
        raise ValueError("The finite root-system data is not initialized")
    finite_weight_space = finite_root_system.weight_space()
    character_ring = WeylCharacterRing(finite_type, style="coroots")
    finite_rho = sum(finite_weight_space.fundamental_weights().values())
    decomposition: dict["AffineWeight", int] = {}

    while remaining:
        dominant_weights = [
            weights_by_key[key]
            for key in remaining
            if weights_by_key[key].finite_part.is_dominant()
        ]
        if not dominant_weights:
            raise ValueError("Weight multiplicities are not finite-Weyl invariant")
        highest_weight = max(
            dominant_weights,
            key=lambda weight: algebra.scalar_product(
                weight.finite_part + finite_rho,
                weight.finite_part + finite_rho,
            ),
        )
        highest_key = affine_weight_key(highest_weight)
        coefficient = remaining[highest_key]
        labels = tuple(highest_weight.finite_dynkin_labels())
        if any(label not in ZZ or label < 0 for label in labels):
            raise ValueError("Irreducible highest weights must be dominant integral")

        decomposition[highest_weight] = coefficient
        for finite_weight, multiplicity in character_ring(labels).weight_multiplicities().items():
            affine_weight = AffineWeight(
                algebra,
                finite_weight_space(finite_weight.to_weight_space()),
                reference_weight.level,
                reference_weight.grade,
            )
            key = affine_weight_key(affine_weight)
            if key not in remaining:
                raise ValueError(
                    "Weight multiplicities do not contain the full finite irreducible "
                    f"character below {highest_weight}"
                )
            remaining[key] -= coefficient * int(multiplicity)
            if remaining[key] < 0:
                raise ArithmeticError(
                    "Finite irreducible subtraction produced a negative multiplicity "
                    f"at {affine_weight}"
                )
            if remaining[key] == 0:
                del remaining[key]
    return decomposition


def compute_bounded_freudenthal_multiplicities(
    tree: BoundedWeightTree,
    *,
    degenerate_multiplicity: Callable[["AffineWeight"], Any],
) -> dict["AffineWeightKey", int]:
    """Compute bounded multiplicities, routing zero denominators to a callback."""
    from .affine_weight import affine_weight_key

    highest_weight = tree.highest_weight
    algebra = highest_weight.algebra
    finite_root_system = algebra._finite_root_system
    if finite_root_system is None:
        raise ValueError("The affine root-system data is not initialized")
    highest_labels = highest_weight.dynkin_labels()
    if any(
        highest_labels[index] not in ZZ or highest_labels[index] < 0
        for index in finite_root_system.index_set()
    ):
        raise NotImplementedError(
            "Finite Weyl reduction requires a highest weight that is dominant "
            "integral for the finite simple roots"
        )
    positive_roots = tuple(
        (
            root,
            _affine_root_from_coordinates(algebra, root.coordinates),
        )
        for root in bounded_affine_positive_roots(algebra, tree.order)
    )
    targets = []
    for weight in tree.weight_by_key.values():
        if classically_dominant_affine_weight(weight) != weight:
            continue
        coordinates = affine_simple_root_coordinates(algebra, highest_weight - weight)
        if any(coordinate not in ZZ or coordinate < 0 for coordinate in coordinates):
            raise ValueError("The bounded weight tree contains a weight outside the highest-weight cone")
        targets.append((sum(coordinates), affine_weight_key(weight), weight, coordinates))
    targets.sort(key=lambda item: (item[0], repr(item[1])))

    highest_key = affine_weight_key(highest_weight)
    multiplicities: dict["AffineWeightKey", int] = {highest_key: 1}
    for _, target_key, target_weight, difference_coordinates in targets:
        if target_key == highest_key:
            continue
        denominator = freudenthal_denominator(highest_weight, target_weight)
        if denominator == 0:
            value = QQ(degenerate_multiplicity(target_weight))
        else:
            numerator = QQ(0)
            for root_record, root_weight in positive_roots:
                max_step = min(
                    difference_coordinate // root_coordinate
                    for difference_coordinate, root_coordinate in zip(
                        difference_coordinates,
                        root_record.coordinates,
                        strict=True,
                    )
                    if root_coordinate > 0
                )
                for step in range(1, int(max_step) + 1):
                    weight_above = target_weight + step * root_weight
                    weight_above_key = affine_weight_key(
                        classically_dominant_affine_weight(weight_above)
                    )
                    multiplicity_above = multiplicities.get(weight_above_key)
                    if multiplicity_above is None:
                        if tree.contains_key(weight_above_key):
                            raise RuntimeError(
                                "A Freudenthal dependency is inside the bounded weight "
                                f"domain but has not been computed: {weight_above}"
                            )
                        continue
                    contribution = (
                        2
                        * root_record.multiplicity
                        * weight_above.scalar_product(root_weight)
                        * multiplicity_above
                    )
                    numerator += contribution
            value = numerator / denominator

        if value not in ZZ or value < 0:
            logger.error(
                "Invalid Freudenthal value target=%s denominator=%s numerator=%s value=%s",
                target_weight,
                denominator,
                numerator if denominator != 0 else None,
                value,
            )
            raise ArithmeticError(
                f"Freudenthal multiplicity must be a nonnegative integer, got {value} "
                f"for {target_weight}"
            )
        multiplicities[target_key] = int(value)
    return multiplicities


class KazhdanLusztigFreudenthalCharacter:
    """Compute a bounded affine character with Freudenthal and lazy KL data."""

    def __init__(
        self,
        algebra: "AffineLieAlgebra",
        *,
        kl_character: Any | None = None,
        orbit_cache_dir: Path | None = None,
    ) -> None:
        from .character import KazhdanLusztigCharacter

        _require_untwisted_affine(algebra)
        self.algebra = algebra
        self.kl_character = kl_character or KazhdanLusztigCharacter(algebra)
        self.orbit_cache_dir = orbit_cache_dir
        if self.kl_character.algebra != algebra:
            raise ValueError("The KL character belongs to a different algebra")
        self._last_stats = HybridCharacterStats()

    def profile_stats(self) -> dict[str, Any]:
        """Return measurements from the most recent multiplicity calculation."""
        return self._last_stats.as_dict()

    def _degenerate_targets(self, tree: BoundedWeightTree) -> tuple["AffineWeight", ...]:
        targets = []
        highest_weight = tree.highest_weight
        for weight in tree.weight_by_key.values():
            if weight == highest_weight:
                continue
            if classically_dominant_affine_weight(weight) != weight:
                continue
            if freudenthal_denominator(highest_weight, weight) == 0:
                targets.append(weight)
        return tuple(targets)

    def multiplicities(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
    ) -> dict[int, dict["AffineWeight", int]]:
        """Return dominant finite-Weyl multiplicities grouped by relative depth."""
        from .affine_weight import affine_weight_key

        if lambda_hat.algebra != self.algebra:
            raise ValueError("The highest weight belongs to a different algebra")

        stats = HybridCharacterStats()
        total_started = time.perf_counter()
        kl_polynomial = self.kl_character.kl_polynomial
        supports_kl_profiling = all(
            hasattr(kl_polynomial, method_name)
            for method_name in ("reset_profile_stats", "set_profiling", "profile_stats")
        )
        if supports_kl_profiling:
            kl_polynomial.reset_profile_stats()
            kl_polynomial.set_profiling(True)

        phase_started = time.perf_counter()
        orbit_arguments: dict[str, Any] = {"order": order}
        if self.orbit_cache_dir is not None:
            orbit_arguments["orbit_cache_dir"] = self.orbit_cache_dir
        orbit = self.kl_character.prepare_bounded_kl_orbit(lambda_hat, **orbit_arguments)
        stats.orbit_seconds = time.perf_counter() - phase_started

        phase_started = time.perf_counter()
        tree = build_bounded_weight_tree(lambda_hat, order=order)
        stats.weight_tree_seconds = time.perf_counter() - phase_started
        degenerate_targets = self._degenerate_targets(tree)
        stats.degenerate_targets = len(degenerate_targets)

        phase_started = time.perf_counter()
        candidate_supports_by_key = {
            affine_weight_key(target): orbit.candidates_above(target)
            for target in degenerate_targets
        }
        candidate_supports = tuple(candidate_supports_by_key.values())
        coordinate_count = self.algebra.rank + 1
        max_coordinates = tuple(
            max(
                (
                    candidate.difference_coordinates[index]
                    for support in candidate_supports
                    for candidate in support
                ),
                default=0,
            )
            for index in range(coordinate_count)
        )
        partition = BoundedAffineVermaPartition(self.algebra, max_coordinates)
        stats.partition_seconds = time.perf_counter() - phase_started

        q_tilde_cache: dict[tuple[int, ...], Any] = {}

        def degenerate_multiplicity(target_weight: "AffineWeight") -> Any:
            return kl_coefficient_at_weight(
                target_weight,
                orbit,
                partition,
                self.kl_character.kl_polynomial,
                q_tilde_cache=q_tilde_cache,
                stats=stats,
                candidates=candidate_supports_by_key[affine_weight_key(target_weight)],
            )

        phase_started = time.perf_counter()
        dominant_multiplicities = compute_bounded_freudenthal_multiplicities(
            tree,
            degenerate_multiplicity=degenerate_multiplicity,
        )
        stats.multiplicity_seconds = time.perf_counter() - phase_started
        stats.freudenthal_targets = len(dominant_multiplicities) - 1
        stats.total_seconds = time.perf_counter() - total_started
        if supports_kl_profiling:
            kl_stats = kl_polynomial.profile_stats()
            stats.kl_q_calls = int(kl_stats["Q_calls"])
            stats.kl_q_cache_hits = int(kl_stats["Q_cache_hits_at_one"]) + int(
                kl_stats["Q_cache_hits_poly"]
            )
            stats.kl_q_tilde_calls = int(kl_stats["Q_tilde_calls"])
            kl_polynomial.set_profiling(False)
        self._last_stats = stats

        return {
            depth: {
                weight: dominant_multiplicities[affine_weight_key(weight)]
                for weight in tree.weights_at_depth(depth)
                if affine_weight_key(weight) in dominant_multiplicities
            }
            for depth in range(tree.order + 1)
        }

    def character(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
    ) -> dict[int, dict["AffineWeight", int]]:
        """Return all bounded weight coefficients grouped by relative depth."""
        dominant_by_depth = self.multiplicities(lambda_hat, order=order)
        tree = build_bounded_weight_tree(lambda_hat, order=order)
        result: dict[int, dict["AffineWeight", int]] = {}
        for depth in range(tree.order + 1):
            dominant_multiplicities = dominant_by_depth[depth]
            result[depth] = {
                weight: dominant_multiplicities[
                    classically_dominant_affine_weight(weight)
                ]
                for weight in tree.weights_at_depth(depth)
            }
        return result

    def finite_irreducible_decomposition(
        self,
        lambda_hat: "AffineWeight",
        *,
        order: int,
    ) -> dict[int, dict["AffineWeight", int]]:
        """Return the finite-dimensional irreducible decomposition at each depth."""
        return {
            depth: decompose_finite_irreducibles(weight_multiplicities)
            for depth, weight_multiplicities in self.character(
                lambda_hat, order=order
            ).items()
        }
