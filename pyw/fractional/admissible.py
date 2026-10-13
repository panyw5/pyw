"""Exact Kac-Wakimoto admissibility for finite and untwisted affine types.

Finite types use all positive coroots; their level argument is metadata.
Affine types use positive *real* coroots, not imaginary roots. The affine
Weyl vector is characterized by its simple-coroot pairings being one.
"""

from fractions import Fraction
from typing import Any

from sage.all import QQ, ZZ, RootSystem, VectorSpace


def _rational(value: Any) -> Any:
    # Sage 10.4 does not coerce Python 3.12 Fraction reliably.
    if isinstance(value, Fraction):
        return QQ((value.numerator, value.denominator))
    if isinstance(value, float):
        return QQ(value.as_integer_ratio())
    return QQ(value)


class AdmissibleWeight:
    """Check regularity and full rational span of integral real coroots.

    For affine input, the supplied level must equal the weight's central
    pairing. Twisted types are explicitly unsupported rather than tested
    using the wrong untwisted real-root formula.
    """

    def __init__(self, cartan_type: Any, weight: Any, level: Any) -> None:
        from pyw.fractional.level import FractionalLevel

        self.cartan_type = tuple(RootSystem(cartan_type).cartan_type())
        self.weight = weight
        if isinstance(level, FractionalLevel):
            if RootSystem(level.cartan_type).cartan_type() != RootSystem(cartan_type).cartan_type():
                raise ValueError("Level and weight Cartan types must agree")
            level = level.level
        try:
            rational = _rational(level)
            self.level = Fraction(int(rational.numerator()), int(rational.denominator()))
        except (TypeError, ValueError) as err:
            raise TypeError("level must be a rational number or FractionalLevel") from err
        self._setup()

    @property
    def root_system(self) -> Any:
        return self._root_system

    @property
    def weight_space(self) -> Any:
        return self._weight_space

    @property
    def weyl_group(self) -> Any:
        """Return the original Weyl-group wrapper, also for finite types."""
        from pyw.core.weyl_group import AffineWeylGroup

        if not hasattr(self, "_weyl_group"):
            self._weyl_group = AffineWeylGroup(self.cartan_type)
        return self._weyl_group

    @property
    def rho_hat(self) -> Any:
        return self._rho_hat

    def _setup(self) -> None:
        self._root_system = RootSystem(self.cartan_type)
        ct = self._root_system.cartan_type()
        if ct.is_affine() and not ct.is_untwisted_affine():
            raise NotImplementedError("Only finite and untwisted affine types are supported")
        self._weight_space = self._root_system.weight_space(QQ, extended=ct.is_affine())
        self._convert_weight_if_needed()
        try:
            self.weight = self._weight_space(self.weight)
        except (TypeError, ValueError) as err:
            raise ValueError("Weight must belong to the specified Cartan type") from err
        if ct.is_affine() and self.weight.level() != _rational(self.level):
            raise ValueError("Weight central level must equal the supplied level")
        self._rho_hat = self._compute_affine_weyl_vector()
        self._cache_coroots()

    def _convert_weight_if_needed(self) -> None:
        from pyw.core.affine_weight import AffineWeight

        if not isinstance(self.weight, AffineWeight):
            return
        aw = self.weight
        if aw.algebra.root_system().cartan_type() != self._root_system.cartan_type():
            raise ValueError("Weight must belong to the specified Cartan type")
        lambdas = self._weight_space.fundamental_weights()
        labels = aw.dynkin_labels()
        self.weight = sum((_rational(value) * lambdas[i] for i, value in labels.items()), self._weight_space.zero())
        # The grade is irrelevant to all real-coroot pairings.

    def _compute_affine_weyl_vector(self) -> Any:
        return self._weight_space.rho()

    def _cache_coroots(self) -> None:
        ct = self._root_system.cartan_type()
        self._finite_root_system = RootSystem(ct.classical()) if ct.is_affine() else self._root_system
        self._coroots_positive = list(self._finite_root_system.coroot_lattice().positive_roots())
        self._coroots_all = self._coroots_positive + [-b for b in self._coroots_positive]

    def _pairing(self, weight: Any, coroot: Any) -> Any:
        return QQ(weight.scalar(coroot))

    def _affine_coroot_strings(self) -> Any:
        """Yield (finite coordinates, initial pairing, central step).

        For a finite root alpha of squared length L (long roots normalized
        to 2), (alpha+n delta)^vee = alpha^vee + n*(2/L)*K.
        Positive strings start at n=0 for alpha>0, at n=1 for alpha<0.
        """
        shifted = self.weight + self._rho_hat
        ct = self._root_system.cartan_type()
        indices = list(self._finite_root_system.cartan_type().index_set())
        finite_roots = list(self._finite_root_system.root_lattice().roots())
        long_norm = max(alpha.norm_squared() for alpha in finite_roots)
        central = [QQ(ct.acheck()[i]) for i in ct.index_set()]
        shifted_level = QQ(shifted.level())
        for b in self._coroots_all:
            coords = [QQ(b.coefficient(i)) for i in indices]
            x = sum(coords[j] * shifted.coefficient(i) for j, i in enumerate(indices))
            factor = QQ(long_norm / b.associated_coroot().norm_squared())
            start = 0 if b.is_positive_root() else 1
            vector = [QQ(0)] + coords
            vector = [v + start * factor * k for v, k in zip(vector, central, strict=True)]
            yield vector, QQ(x + start * factor * shifted_level), factor * shifted_level, [factor * k for k in central]

    @staticmethod
    def _integral_residue(x: Any, step: Any) -> int | None:
        """First n>=0 such that x+n*step is integral, or None.

        Rational steps have period denominator(step); solve the congruence
        directly, avoiding a potentially enormous enumeration.
        """
        if step == 0:
            return 0 if x in ZZ else None
        period = step.denominator()
        scaled = x * period
        if scaled not in ZZ:
            return None
        return int((-ZZ(scaled) * step.numerator().inverse_mod(period)) % period) if period != 1 else 0

    def _check_condition_1a(self) -> bool:
        if not self._root_system.cartan_type().is_affine():
            shifted = self.weight + self._rho_hat
            return all(not (self._pairing(shifted, b) in ZZ and self._pairing(shifted, b) <= 0) for b in self._coroots_positive)
        for _, x, step, _ in self._affine_coroot_strings():
            n = self._integral_residue(x, step)
            if n is not None and (step < 0 or x + n * step <= 0):
                return False
        return True

    def _check_condition_1b(self) -> bool:
        if not self._root_system.cartan_type().is_affine():
            shifted = self.weight + self._rho_hat
            integral = [b for b in self._coroots_positive if self._pairing(shifted, b) in ZZ]
            return self._compute_q_span(integral) == self._compute_q_span(self._coroots_all)
        space = VectorSpace(QQ, len(self._root_system.cartan_type().index_set()))
        vectors = []
        for base, x, step, direction in self._affine_coroot_strings():
            n = self._integral_residue(x, step)
            if n is not None:
                vectors.append(space([v + n * d for v, d in zip(base, direction, strict=True)]))
                # Two successive integral roots give the central direction.
                vectors.append(space(direction))
        return space.span(vectors).dimension() == space.dimension()

    def _compute_q_span(self, coroots: list[Any]) -> Any:
        lattice = self._finite_root_system.coroot_lattice()
        space = VectorSpace(QQ, len(lattice.index_set()))
        return space.span([space([b.coefficient(i) for i in lattice.index_set()]) for b in coroots])

    def _spans_equal(self, span1: Any, span2: Any) -> bool:
        return bool(span1 == span2)

    def is_admissible(self) -> bool:
        return bool(self._check_condition_1a() and self._check_condition_1b())

    def __repr__(self) -> str:
        return f"AdmissibleWeight({self.cartan_type}, admissible={self.is_admissible()})"

    def __str__(self) -> str:
        status = "admissible" if self.is_admissible() else "not admissible"
        return f"AdmissibleWeight({self.cartan_type}, status={status})"
