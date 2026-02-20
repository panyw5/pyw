"""
Formal Character Module for Affine Lie Algebras

This module provides classes for representing and manipulating formal characters
of modules over affine Lie algebras, including:

- Truncated formal power series in q = e^{-δ}
- Weyl-Kac denominator formula
- Character arithmetic (addition, multiplication, scalar)

The character of a module M is represented as:
    ch(M) = Σ_{n} a_n q^n

where a_n is the multiplicity (or weight space dimension) at grade n.

References:
    - Kac, V. G. "Infinite-dimensional Lie algebras" (3rd ed.), Chapter 10
    - Di Francesco et al., "Conformal Field Theory", Chapter 14
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional, Tuple, Union

from sage.all import QQ, Integer, prod

if TYPE_CHECKING:
    from .affine_lie_algebra import AffineLieAlgebra
    from .affine_weight import AffineWeight


@dataclass
class FormalCharacter:
    """
    A formal character as a truncated power series in q = e^{-δ}.

    The character is represented as a dictionary mapping grades (powers of q)
    to coefficients (multiplicities or weight space contributions).

    Parameters
    ----------
    coefficients : Dict[int, Any]
        Mapping from grade n to coefficient a_n
        Represents: ch = Σ a_n q^n
    max_grade : int
        Maximum grade computed (truncation bound)
    algebra : AffineLieAlgebra, optional
        The affine Lie algebra context

    Examples
    --------
    >>> # Character 1 + 2q + 3q^2
    >>> ch = FormalCharacter({0: 1, 1: 2, 2: 3}, max_grade=2)
    >>> ch[0]  # Coefficient of q^0
    1
    >>> ch[1]  # Coefficient of q^1
    2

    Notes
    -----
    The grade n corresponds to the L₀ eigenvalue in the Di Francesco notation:
    a state with affine weight (λ; k; n) contributes to the q^n term.
    """

    coefficients: Dict[int, Any] = field(default_factory=dict)
    max_grade: int = 10
    algebra: Optional["AffineLieAlgebra"] = None

    def __post_init__(self):
        """Normalize coefficients."""
        # Remove zero coefficients
        self.coefficients = {
            n: c for n, c in self.coefficients.items() if c != 0 and n <= self.max_grade
        }

    def __getitem__(self, grade: int) -> Any:
        """Get coefficient at given grade."""
        return self.coefficients.get(grade, 0)

    def __setitem__(self, grade: int, value: Any) -> None:
        """Set coefficient at given grade."""
        if value == 0:
            self.coefficients.pop(grade, None)
        elif grade <= self.max_grade:
            self.coefficients[grade] = value

    def __iter__(self) -> Iterator[Tuple[int, Any]]:
        """Iterate over (grade, coefficient) pairs."""
        for n in sorted(self.coefficients.keys()):
            yield n, self.coefficients[n]

    def __len__(self) -> int:
        """Number of non-zero terms."""
        return len(self.coefficients)

    # =========================================================================
    # Arithmetic Operations
    # =========================================================================

    def __add__(self, other: "FormalCharacter") -> "FormalCharacter":
        """Add two characters."""
        if not isinstance(other, FormalCharacter):
            return NotImplemented

        max_grade = max(self.max_grade, other.max_grade)
        result = {}

        all_grades = set(self.coefficients.keys()) | set(other.coefficients.keys())
        for n in all_grades:
            if n <= max_grade:
                coeff = self[n] + other[n]
                if coeff != 0:
                    result[n] = coeff

        return FormalCharacter(result, max_grade, self.algebra or other.algebra)

    def __sub__(self, other: "FormalCharacter") -> "FormalCharacter":
        """Subtract two characters."""
        if not isinstance(other, FormalCharacter):
            return NotImplemented

        max_grade = max(self.max_grade, other.max_grade)
        result = {}

        all_grades = set(self.coefficients.keys()) | set(other.coefficients.keys())
        for n in all_grades:
            if n <= max_grade:
                coeff = self[n] - other[n]
                if coeff != 0:
                    result[n] = coeff

        return FormalCharacter(result, max_grade, self.algebra or other.algebra)

    def __mul__(self, other: Union["FormalCharacter", int, Any]) -> "FormalCharacter":
        """Multiply character by scalar or another character."""
        if isinstance(other, (int, Integer)) or hasattr(other, "__rmul__"):
            # Scalar multiplication
            result = {n: c * other for n, c in self.coefficients.items() if c * other != 0}
            return FormalCharacter(result, self.max_grade, self.algebra)

        if isinstance(other, FormalCharacter):
            # Character multiplication (convolution)
            max_grade = min(self.max_grade, other.max_grade)
            result = {}

            for n1, c1 in self.coefficients.items():
                for n2, c2 in other.coefficients.items():
                    n = n1 + n2
                    if n <= max_grade:
                        result[n] = result.get(n, 0) + c1 * c2

            # Remove zeros
            result = {n: c for n, c in result.items() if c != 0}
            return FormalCharacter(result, max_grade, self.algebra or other.algebra)

        return NotImplemented

    def __rmul__(self, scalar: Any) -> "FormalCharacter":
        """Right scalar multiplication."""
        return self.__mul__(scalar)

    def __neg__(self) -> "FormalCharacter":
        """Negate character."""
        return self * (-1)

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def leading_grade(self) -> Optional[int]:
        """The smallest grade with non-zero coefficient."""
        if not self.coefficients:
            return None
        return min(self.coefficients.keys())

    @property
    def leading_coefficient(self) -> Any:
        """The coefficient of the leading term."""
        if self.leading_grade is None:
            return 0
        return self.coefficients[self.leading_grade]

    def truncate(self, new_max: int) -> "FormalCharacter":
        """Return a truncated copy."""
        result = {n: c for n, c in self.coefficients.items() if n <= new_max}
        return FormalCharacter(result, new_max, self.algebra)

    def shift(self, delta: int) -> "FormalCharacter":
        """Shift all grades by delta (multiply by q^delta)."""
        result = {n + delta: c for n, c in self.coefficients.items()}
        return FormalCharacter(result, self.max_grade + delta, self.algebra)

    # =========================================================================
    # Representation
    # =========================================================================

    def __repr__(self) -> str:
        if not self.coefficients:
            return "FormalCharacter(0)"
        terms = []
        for n in sorted(self.coefficients.keys())[:5]:
            c = self.coefficients[n]
            if n == 0:
                terms.append(str(c))
            elif n == 1:
                terms.append(f"{c}*q")
            else:
                terms.append(f"{c}*q^{n}")
        result = " + ".join(terms)
        if len(self.coefficients) > 5:
            result += " + ..."
        return f"FormalCharacter({result})"

    def _repr_latex_(self) -> str:
        """LaTeX representation for Jupyter."""
        if not self.coefficients:
            return "$0$"
        terms = []
        for n in sorted(self.coefficients.keys())[:10]:
            c = self.coefficients[n]
            if c == 1:
                c_str = ""
            elif c == -1:
                c_str = "-"
            else:
                c_str = str(c)

            if n == 0:
                terms.append(str(c))
            elif n == 1:
                terms.append(f"{c_str}q")
            else:
                terms.append(f"{c_str}q^{{{n}}}")

        result = " + ".join(terms).replace("+ -", "- ")
        if len(self.coefficients) > 10:
            result += " + \\cdots"
        return f"${result}$"

    def to_dict(self) -> Dict[int, Any]:
        """Return coefficients as a dictionary."""
        return dict(self.coefficients)

    def to_list(self, up_to: Optional[int] = None) -> list:
        """Return coefficients as a list [a_0, a_1, ..., a_n]."""
        if up_to is None:
            up_to = self.max_grade
        return [self[n] for n in range(up_to + 1)]


class WeylKacDenominator:
    """
    Weyl-Kac denominator formula for affine Lie algebras.

    The denominator is:
        Π_{α > 0} (1 - e^{-α})^{mult(α)}

    which by the Weyl-Kac formula equals:
        Σ_{w ∈ W} (-1)^{ℓ(w)} e^{w(ρ̂) - ρ̂}

    For character computations, we need the inverse:
        1 / denominator = Σ_n d_n q^n

    Parameters
    ----------
    algebra : AffineLieAlgebra
        The affine Lie algebra

    Examples
    --------
    >>> from pyw.core import AffineLieAlgebra
    >>> ala = AffineLieAlgebra(['A', 2, 1])
    >>> denom = WeylKacDenominator(ala)
    >>> inv = denom.inverse(max_grade=5)
    >>> inv[0]  # Leading coefficient
    1

    References
    ----------
    Di Francesco et al., Eq. (14.150)
    """

    def __init__(self, algebra: "AffineLieAlgebra") -> None:
        self._algebra = algebra
        self._inverse_cache: Dict[int, FormalCharacter] = {}

    @property
    def algebra(self) -> "AffineLieAlgebra":
        """The affine Lie algebra."""
        return self._algebra

    def inverse(self, max_grade: int = 10) -> FormalCharacter:
        """
        Compute the inverse denominator as a q-series.

        The inverse denominator appears in character formulas:
            ch(M(λ)) = e^λ / denominator

        Parameters
        ----------
        max_grade : int
            Maximum grade to compute

        Returns
        -------
        FormalCharacter
            The inverse denominator truncated at max_grade
        """
        # Check cache
        if max_grade in self._inverse_cache:
            return self._inverse_cache[max_grade]

        # Compute using product formula
        # 1/denominator = Π_{α > 0} (1 - q^{n_α})^{-mult(α)}
        #               = Π_{α > 0} Σ_{k≥0} p(k) q^{k·n_α}
        # where p(k) is the partition function contribution

        result = self._compute_inverse_product(max_grade)
        self._inverse_cache[max_grade] = result
        return result

    def _compute_inverse_product(self, max_grade: int) -> FormalCharacter:
        """
        Compute inverse denominator using the product formula.

        For affine algebras, the denominator factors as:
            denominator = η(q)^r · (finite Weyl denominator terms)

        where η(q) = q^{1/24} Π_{n≥1} (1-q^n) is the Dedekind eta function
        and r is the rank.
        """
        rank = self._algebra.rank if hasattr(self._algebra, "rank") else 1

        # Start with 1
        coeffs = {0: QQ(1)}

        # For each positive root contribution
        # Simplified: use partition function approximation
        # 1/(1-q^n) = 1 + q^n + q^{2n} + ...

        for n in range(1, max_grade + 1):
            # Contribution from imaginary roots (multiplicity = rank)
            # 1/(1-q^n)^rank
            new_coeffs = {}
            for grade, coeff in coeffs.items():
                # Add contributions from (1-q^n)^{-rank}
                for k in range((max_grade - grade) // n + 1):
                    new_grade = grade + k * n
                    if new_grade <= max_grade:
                        # Binomial coefficient for (1-x)^{-rank}
                        binom = self._negative_binomial(rank, k)
                        new_coeffs[new_grade] = new_coeffs.get(new_grade, 0) + coeff * binom

            coeffs = new_coeffs

        return FormalCharacter(coeffs, max_grade, self._algebra)

    def _negative_binomial(self, r: int, k: int) -> Any:
        """
        Compute binomial coefficient C(-r, k) = C(r+k-1, k) * (-1)^k.

        For (1-x)^{-r} = Σ_k C(r+k-1, k) x^k
        """
        from sage.all import binomial

        return binomial(r + k - 1, k)


class VermaCharacter:
    """
    Character of a Verma module M(λ).

    The Verma module character is:
        ch(M(λ)) = e^λ / Π_{α > 0} (1 - e^{-α})^{mult(α)}
                 = e^λ · (inverse denominator)

    In q-expansion form (q = e^{-δ}):
        ch(M(λ)) = q^{-n_λ} · Σ_k c_k q^k

    where n_λ is the grade of λ.

    Parameters
    ----------
    algebra : AffineLieAlgebra
        The affine Lie algebra
    weight : AffineWeight
        The highest weight λ

    Examples
    --------
    >>> from pyw.core import AffineLieAlgebra, AffineWeight
    >>> ala = AffineLieAlgebra(['A', 2, 1])
    >>> Lambda = ala.fundamental_weights()
    >>> verma = VermaCharacter(ala, Lambda[1])
    >>> ch = verma.character(max_grade=5)
    """

    def __init__(self, algebra: "AffineLieAlgebra", weight: "AffineWeight") -> None:
        self._algebra = algebra
        self._weight = weight
        self._denominator = WeylKacDenominator(algebra)

    @property
    def algebra(self) -> "AffineLieAlgebra":
        """The affine Lie algebra."""
        return self._algebra

    @property
    def weight(self) -> "AffineWeight":
        """The highest weight."""
        return self._weight

    def character(self, max_grade: int = 10) -> FormalCharacter:
        """
        Compute the Verma module character.

        Parameters
        ----------
        max_grade : int
            Maximum grade to compute

        Returns
        -------
        FormalCharacter
            The character ch(M(λ)) truncated at max_grade
        """
        # Get the grade (n-value) of the weight
        weight_grade = int(self._weight.grade) if hasattr(self._weight, "grade") else 0

        # Get inverse denominator
        inv_denom = self._denominator.inverse(max_grade + abs(weight_grade))

        # Shift by weight grade: ch(M(λ)) = q^{-n_λ} · inv_denom
        # Note: grade n corresponds to q^n, so e^λ with grade n_λ gives q^{n_λ}
        result = inv_denom.shift(-weight_grade)

        # Truncate to requested max_grade
        return result.truncate(max_grade)


def character_from_weight(
    algebra: "AffineLieAlgebra",
    weight: "AffineWeight",
    max_grade: int = 10,
) -> FormalCharacter:
    """
    Convenience function to compute Verma character.

    Parameters
    ----------
    algebra : AffineLieAlgebra
        The affine Lie algebra
    weight : AffineWeight
        The highest weight
    max_grade : int
        Maximum grade

    Returns
    -------
    FormalCharacter
        The Verma module character
    """
    verma = VermaCharacter(algebra, weight)
    return verma.character(max_grade)
