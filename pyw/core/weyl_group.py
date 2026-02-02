"""
Affine Weyl Group Wrapper Module

This module provides a clean wrapper around SageMath's Weyl group implementation
for affine types, including support for finite, affine, and extended affine Weyl groups.

Key Features:
    - Finite Weyl groups (e.g., type A2)
    - Affine Weyl groups (e.g., type A2^(1))
    - Extended affine Weyl groups
    - Dot action: w.λ = w(λ + ρ̂) - ρ̂
    - Translation operators: t_α

References:
    - Kac, V. G., "Infinite-dimensional Lie algebras" (3rd ed.)
    - Kac, V. G., Wakimoto, M. "On rationality of W-algebras"
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Iterable, Iterator, Mapping, Sequence, Union

from sage.all import QQ, RootSystem, WeylGroup


class AffineWeylGroup:
    """
    Wrapper class for SageMath's Weyl group with special support for affine types.

    This class provides a convenient interface for working with finite, affine,
    and extended affine Weyl groups, including the dot action and translation
    operators essential for admissible weight calculations.

    Parameters
    ----------
    cartan_type : Union[str, tuple]
        Cartan type specification. Examples:
        - 'A2' for finite type A2
        - ['A', 2, 1] for affine type A2^(1)
        - ('A', 2, 1) for affine type A2^(1)
    extended : bool, optional (default=False)
        If True, use the extended affine Weyl group (Ŵ = W ⋉ t_Q*).
        If False, use the standard affine Weyl group (Ŵ = W ⋉ t_Q∨).
        For finite types, this parameter has no effect.

    Attributes
    ----------
    cartan_type : tuple
        The Cartan type specification.
    extended : bool
        Whether the extended affine Weyl group is used.

    Examples
    --------
    >>> from pyw.core.weyl_group import AffineWeylGroup
    >>>
    >>> # Finite Weyl group
    >>> W_finite = AffineWeylGroup('A2')
    >>> print(W_finite.weyl_group)
    Weyl Group of type ['A', 2]
    >>>
    >>> # Affine Weyl group
    >>> W_affine = AffineWeylGroup(['A', 2, 1])
    >>> print(W_affine.weyl_group)
    Weyl Group of type ['A', 2, 1]
    >>>
    >>> # Extended affine Weyl group
    >>> W_extended = AffineWeylGroup(['A', 2, 1], extended=True)
    >>>
    >>> # Dot action
    >>> w = W_affine.weyl_group.simple_reflection(0)
    >>> weight = W_affine.weyl_group.domain().fundamental_weights()[0]
    >>> rho = W_affine.weyl_group.domain().rho()
    >>> result = W_affine.dot_action(w, weight, rho)
    """

    def __init__(self, cartan_type: Union[str, tuple, list], extended: bool = False) -> None:
        """
        Initialize the AffineWeylGroup wrapper.

        Parameters
        ----------
        cartan_type : Union[str, tuple, list]
            Cartan type specification.
        extended : bool, optional (default=False)
            Whether to use extended affine Weyl group.

        Raises
        ------
        ValueError
            If the cartan_type is not recognized by SageMath.
        """
        self.cartan_type = cartan_type
        self.extended = extended
        self._weyl_group: Any = None
        self._simple_reflections: dict[int, Any] = {}
        self._setup()

    def _setup(self) -> None:
        """
        Initialize the SageMath Weyl group and related structures.

        This method creates the underlying SageMath WeylGroup object
        and caches the simple reflections for efficient access.
        """
        # Create the Weyl group using SageMath
        self._weyl_group = WeylGroup(self.cartan_type)

        # Cache simple reflections
        for i in self._weyl_group.index_set():
            self._simple_reflections[i] = self._weyl_group.simple_reflection(i)

    @property
    def weyl_group(self) -> Any:
        """
        Get the underlying SageMath WeylGroup object.

        Returns
        -------
        WeylGroup
            The SageMath WeylGroup object for the specified Cartan type.

        Examples
        --------
        >>> W = AffineWeylGroup(['A', 2, 1])
        >>> wg = W.weyl_group
        >>> print(wg)
        Weyl Group of type ['A', 2, 1]
        """
        return self._weyl_group

    @property
    def simple_reflections(self) -> dict[int, Any]:
        """
        Get the simple reflections of the Weyl group.

        Returns
        -------
        dict[int, Any]
            Dictionary mapping node indices to simple reflections.
            For affine types, the mapping includes the affine node (index -1 or len-1).

        Examples
        --------
        >>> W = AffineWeylGroup(['A', 2, 1])
        >>> s0 = W.simple_reflections[0]
        >>> s1 = W.simple_reflections[1]
        >>> s2 = W.simple_reflections[2]
        """
        return self._simple_reflections

    def dot_action(self, w: Any, weight: Any, rho: Any) -> Any:
        """
        Compute the dot action of a Weyl group element on a weight.

        The dot action is defined as:
            w.λ = w(λ + ρ̂) - ρ̂

        This action is essential for admissible weight calculations and
        appears in the Kac-Wakimoto theory of fractional levels.

        Parameters
        ----------
        w : Any
            A Weyl group element (from `self.weyl_group`).
        weight : Any
            A weight in the weight lattice or weight space.
        rho : Any
            The Weyl vector (ρ̂). For affine types, this is the affine Weyl vector.

        Returns
        -------
        Any
            The result of the dot action w.λ.

        Examples
        --------
        >>> from pyw.core.weyl_group import AffineWeylGroup
        >>> from sage.all import RootSystem
        >>>
        >>> W = AffineWeylGroup(['A', 2, 1])
        >>> R = RootSystem(['A', 2, 1])
        >>> P = R.weight_space(extended=True)
        >>> Lambda = P.fundamental_weights()
        >>>
        >>> # Get a Weyl group element and weight
        >>> s0 = W.simple_reflections[0]
        >>> lambda_weight = Lambda[0]
        >>> rho_hat = P.rho()
        >>>
        >>> # Compute dot action
        >>> result = W.dot_action(s0, lambda_weight, rho_hat)
        >>> print(result)
        """
        # Compute w(λ + ρ)
        shifted_weight = weight + rho
        w_shifted = w.action(shifted_weight)

        # Compute w.λ = w(λ + ρ) - ρ
        result = w_shifted - rho

        return result

    def translation_operator(self, alpha: Any) -> Any:
        """
        Create a translation operator t_α in the extended affine Weyl group.

        Translation operators are elements of the extended affine Weyl group
        that act by translation on the weight lattice. For the affine Weyl
        group Ŵ = W ⋉ t_Q∨, translations are by coroots. For the extended
        affine Weyl group Ŵ = W ⋉ t_Q*, translations are by weights.

        Parameters
        ----------
        alpha : Any
            A root or weight defining the translation direction.
            For standard affine groups: α ∈ Q∨ (coroot lattice)
            For extended affine groups: α ∈ Q* (weight lattice)

        Returns
        -------
        Any
            The translation operator t_α in the Weyl group.

        Raises
        ------
        NotImplementedError
            If translation operators are not directly available for the
            specified Cartan type.

        Examples
        --------
        >>> from pyw.core.weyl_group import AffineWeylGroup
        >>> from sage.all import RootSystem
        >>>
        >>> W = AffineWeylGroup(['A', 2, 1])
        >>> R = RootSystem(['A', 2, 1])
        >>> simple_roots = R.root_lattice().simple_roots()
        >>>
        >>> # Create translation operator (note: implementation depends on SageMath)
        >>> try:
        ...     t_alpha = W.translation_operator(simple_roots[0])
        ...     print(f"Translation: {t_alpha}")
        ... except NotImplementedError:
        ...     print("Translation not directly available")
        """
        # For affine Weyl groups, translation operators are part of the semidirect product
        # The exact implementation depends on SageMath's internal representation

        # Check if we can construct translation using the weight/weight lattice
        if hasattr(self._weyl_group, "translation"):
            # SageMath may provide a translation method
            return self._weyl_group.translation(alpha)
        else:
            # Alternative: Construct translation using the algebraic group structure
            # This requires more advanced SageMath operations
            raise NotImplementedError(
                f"Translation operators for Cartan type {self.cartan_type} "
                "require specialized construction. "
                "Use the underlying SageMath WeylGroup directly for advanced operations."
            )

    def __repr__(self) -> str:
        """
        Return a string representation of the AffineWeylGroup.

        Returns
        -------
        str
            String representation showing Cartan type and extended status.
        """
        extended_str = "extended" if self.extended else "standard"
        return f"AffineWeylGroup({self.cartan_type}, {extended_str})"

    def __str__(self) -> str:
        """
        Return a user-friendly string representation.

        Returns
        -------
        str
            User-friendly string showing the underlying Weyl group.
        """
        return str(self._weyl_group)


def _reduced_word_str(word: Sequence[int]) -> str:
    if not word:
        return "1"
    return " * ".join(f"s_{i}" for i in word)


def _sage_element_length(w: Any) -> int:
    if hasattr(w, "length"):
        try:
            return int(w.length())
        except Exception:
            pass
    if hasattr(w, "reduced_word"):
        try:
            return len(w.reduced_word())
        except Exception:
            pass
    return 0


class FiniteWeylGroup:
    """Finite Weyl group wrapper with elements acting on AffineWeight.

    Notes
    -----
    This is intentionally minimal and only implements the surface used by:
    - pyw/core/affine_lie_algebra.py (finite_weyl_group)
    - debug notebooks/scripts in this repo.
    """

    def __init__(self, algebra: Any) -> None:
        self._algebra = algebra
        finite_rs = (
            algebra._finite_root_system
            if getattr(algebra, "is_affine", False)
            else algebra.root_system()
        )
        # Matrix Weyl group acting on finite weight space.
        self._weight_space = finite_rs.weight_space(QQ)
        self._W = self._weight_space.weyl_group()

    @property
    def algebra(self) -> Any:
        return self._algebra

    def __len__(self) -> int:
        return len(self._W)

    def __iter__(self) -> Iterator[FiniteWeylGroupElement]:
        elts = list(self._W)
        elts_sorted = sorted(elts, key=lambda w: (_sage_element_length(w), tuple(w.reduced_word())))
        for w in elts_sorted:
            yield FiniteWeylGroupElement(self, w)

    def identity(self) -> "FiniteWeylGroupElement":
        return FiniteWeylGroupElement(self, self._W.one())

    def longest_element(self) -> "FiniteWeylGroupElement":
        # Sage uses long_element() for the longest element.
        return FiniteWeylGroupElement(self, self._W.long_element())

    def simple_reflection(self, i: int) -> "FiniteWeylGroupElement":
        return FiniteWeylGroupElement(self, self._W.simple_reflection(i))


@dataclass(frozen=True)
class FiniteWeylGroupElement:
    _group: FiniteWeylGroup
    _w: Any

    def reduced_word(self) -> list[int]:
        return list(self._w.reduced_word())

    def __str__(self) -> str:
        return _reduced_word_str(self.reduced_word())

    def __repr__(self) -> str:
        return f"FiniteWeylGroupElement({str(self)})"

    def __mul__(self, other: "FiniteWeylGroupElement") -> "FiniteWeylGroupElement":
        if self._group is not other._group:
            raise TypeError("Cannot multiply elements from different FiniteWeylGroup instances")
        return FiniteWeylGroupElement(self._group, self._w * other._w)

    def action(self, weight: Any) -> Any:
        # Avoid import cycles.
        from .affine_weight import AffineWeight

        if isinstance(weight, AffineWeight):
            new_finite = self._w.action(weight.finite_part)
            return AffineWeight(weight.algebra, new_finite, weight.level, weight.grade)

        return self._w.action(weight)

    __call__ = action

    @property
    def sage_element(self) -> Any:
        return self._w


class ExtendedAffineWeylGroup:
    """Extended affine Weyl group implemented as a semidirect product W ⋉ Q^∨.

    Representation convention
    -------------------------
    An element is represented as (w, beta) meaning **w * t_beta**.
    Action on an AffineWeight is:
        (w, beta).action(λ) = w.action( t_beta(λ) )
    Multiplication is:
        (w1, b1) * (w2, b2) = (w1*w2, w2^{-1}(b1) + b2)
    """

    def __init__(self, algebra: Any):
        if not getattr(algebra, "is_affine", False):
            raise ValueError("ExtendedAffineWeylGroup requires an affine algebra")

        self._algebra = algebra
        finite_rs = algebra._finite_root_system

        self._finite_weight_space = finite_rs.weight_space(QQ)
        self._finite_root_lattice = finite_rs.root_lattice()
        self._finite_coroot_space = finite_rs.coroot_space(QQ)

        self._W_weight = self._finite_weight_space.weyl_group()
        self._W_root = self._finite_root_lattice.weyl_group()
        self._W_coroot = self._finite_coroot_space.weyl_group()

        self._zero_beta = self._finite_coroot_space.zero()

        self._theta_coroot = self._compute_theta_coroot()
        self._s_theta_weight = self._compute_s_theta_on_weight_space()

    def _translate_affine_weight(self, beta: Any, weight: Any) -> Any:
        """Translate an AffineWeight by a coroot beta.

        We intentionally implement this here (instead of delegating to
        AffineLieAlgebra.translation) because that method uses the library's
        Kac–Moody scalar_product normalization, while the translation formula
        needs the *weight–coroot pairing* <λ, β> and the coroot norm |β|^2.

        The implementation below computes both via Sage's ambient space model,
        which satisfies <Λ_i, α_j^∨> = δ_{ij}.
        """

        from .affine_weight import AffineWeight

        if not isinstance(weight, AffineWeight):
            raise TypeError("Translation action expects an AffineWeight")

        lam = weight.finite_part
        k = weight.level
        n = weight.grade

        # Fast path via Sage ambient space.
        try:
            lam_amb = lam.to_ambient()
            beta_amb = beta.to_ambient()
            v_amb = lam_amb + k * beta_amb

            parent = lam.parent()
            if not hasattr(parent, "fundamental_weights"):
                raise AttributeError("finite_part parent has no fundamental_weights")
            Lambda = parent.fundamental_weights()

            A = beta_amb.parent()
            simple_coroots = A.simple_coroots()

            new_lam = parent.zero()
            for i, Lambda_i in Lambda.items():
                # Defensive: some parents may include index 0 in affine contexts.
                if int(i) == 0:
                    continue
                coeff = v_amb.inner_product(simple_coroots[int(i)])
                if coeff:
                    new_lam += coeff * Lambda_i

            pairing = lam_amb.inner_product(beta_amb)
            if hasattr(beta, "norm_squared"):
                beta_norm_sq = beta.norm_squared()
            else:
                beta_norm_sq = beta_amb.inner_product(beta_amb)

            new_n = n - pairing - k * QQ(beta_norm_sq) / 2
            return AffineWeight(weight.algebra, new_lam, k, new_n)
        except Exception:
            # Fallback to existing implementation if ambient conversions fail.
            return self._algebra.translation(beta, weight)

    @property
    def algebra(self) -> Any:
        return self._algebra

    def __repr__(self) -> str:
        return f"ExtendedAffineWeylGroup({self._algebra._cartan_type})"

    def identity(self) -> "ExtendedAffineWeylGroupElement":
        return ExtendedAffineWeylGroupElement(self, self._W_weight.one(), self._zero_beta)

    def theta_coroot(self) -> Any:
        return self._theta_coroot

    def translation(self, beta: Any) -> "ExtendedAffineWeylGroupElement":
        beta_cs = self._coerce_to_coroot_space(beta)
        return ExtendedAffineWeylGroupElement(self, self._W_weight.one(), beta_cs)

    def simple_reflection(self, i: int) -> "ExtendedAffineWeylGroupElement":
        if i == 0:
            # Convention used by our debug notebook/scripts:
            # s0 = s_theta * t_{-theta^vee}
            return ExtendedAffineWeylGroupElement(self, self._s_theta_weight, -self._theta_coroot)

        return ExtendedAffineWeylGroupElement(
            self, self._W_weight.simple_reflection(i), self._zero_beta
        )

    def elements_as_semi_direct_product(
        self, *, translation_bounds: Mapping[int, tuple[int, int]]
    ) -> list["ExtendedAffineWeylGroupElement"]:
        """Enumerate a bounded subset of W ⋉ Q^∨.

        Parameters
        ----------
        translation_bounds:
            Dict mapping simple coroot indices (1..rank) to (min,max) coefficients.
            Unspecified indices default to (0,0).

        Returns
        -------
        list[ExtendedAffineWeylGroupElement]
            Deterministic ordering: first all finite elements with zero translation,
            then increasing translation (lexicographic coefficients) combined with
            finite elements sorted by (length, reduced_word).
        """

        idxs = list(self._finite_coroot_space.index_set())
        ranges: list[range] = []
        for i in idxs:
            lo, hi = translation_bounds.get(int(i), (0, 0))
            ranges.append(range(int(lo), int(hi) + 1))

        finite_elts = list(self._W_weight)
        finite_elts_sorted = sorted(
            finite_elts, key=lambda w: (_sage_element_length(w), tuple(w.reduced_word()))
        )

        result: list[ExtendedAffineWeylGroupElement] = []
        # Zero translation first.
        for w in finite_elts_sorted:
            result.append(ExtendedAffineWeylGroupElement(self, w, self._zero_beta))

        # Non-zero translations.
        basis = self._finite_coroot_space.simple_roots()  # alphacheck[i]
        for coeffs in product(*ranges):
            if all(c == 0 for c in coeffs):
                continue
            beta = self._zero_beta
            for i, c in zip(idxs, coeffs, strict=True):
                if c:
                    beta += int(c) * basis[int(i)]
            for w in finite_elts_sorted:
                result.append(ExtendedAffineWeylGroupElement(self, w, beta))

        return result

    def _coerce_to_coroot_space(self, beta: Any) -> Any:
        if beta in (0, None):
            return self._zero_beta
        if getattr(beta, "parent", lambda: None)() is self._finite_coroot_space:
            return beta
        try:
            return self._finite_coroot_space(beta)
        except Exception:
            # Best-effort reconstruction from monomial coefficients.
            if hasattr(beta, "monomial_coefficients"):
                coeffs = beta.monomial_coefficients()
                basis = self._finite_coroot_space.simple_roots()
                out = self._zero_beta
                for i, c in coeffs.items():
                    out += c * basis[int(i)]
                return out
            raise

    def _compute_theta_coroot(self) -> Any:
        theta = self._finite_root_lattice.highest_root()
        theta_vee = theta.associated_coroot()
        return self._finite_coroot_space(theta_vee)

    def _compute_s_theta_on_weight_space(self) -> Any:
        # Compute reflection s_theta in the finite Weyl group by conjugating
        # a simple reflection: if u(α_j) = θ, then s_θ = u s_j u^{-1}.
        theta = self._finite_root_lattice.highest_root()
        simple_roots = self._finite_root_lattice.simple_roots()

        candidates = list(self._W_root)
        for j, alpha_j in simple_roots.items():
            for u in candidates:
                if u.action(alpha_j) == theta:
                    u_word = list(u.reduced_word())
                    u_weight = self._W_weight.from_reduced_word(u_word)
                    s_j_weight = self._W_weight.simple_reflection(int(j))
                    return u_weight * s_j_weight * u_weight.inverse()

        raise RuntimeError("Failed to construct s_theta in finite Weyl group")


@dataclass(frozen=True)
class ExtendedAffineWeylGroupElement:
    _group: ExtendedAffineWeylGroup
    _w: Any  # finite Weyl element acting on finite weight space
    _beta: Any  # coroot-space element

    @property
    def finite_part(self) -> FiniteWeylGroupElement:
        return FiniteWeylGroupElement(FiniteWeylGroup(self._group.algebra), self._w)

    @property
    def translation_vector(self) -> Any:
        return self._beta

    def reduced_word(self) -> list[int]:
        return list(self._w.reduced_word())

    def __str__(self) -> str:
        w_str = _reduced_word_str(self.reduced_word())
        if self._beta == self._group._zero_beta:
            return w_str
        if not self.reduced_word():
            return f"t_{{{self._beta}}}"
        return f"{w_str} * t_{{{self._beta}}}"

    def __repr__(self) -> str:
        return f"ExtendedAffineWeylGroupElement({str(self)})"

    def __mul__(self, other: "ExtendedAffineWeylGroupElement") -> "ExtendedAffineWeylGroupElement":
        if self._group is not other._group:
            raise TypeError(
                "Cannot multiply elements from different ExtendedAffineWeylGroup instances"
            )

        # (w1, b1) * (w2, b2) = (w1*w2, w2^{-1}(b1) + b2)
        w_new = self._w * other._w
        word2 = list(other._w.reduced_word())
        w2_coroot = self._group._W_coroot.from_reduced_word(word2)
        b_new = w2_coroot.inverse().action(self._beta) + other._beta
        return ExtendedAffineWeylGroupElement(self._group, w_new, b_new)

    def action(self, weight: Any) -> Any:
        from .affine_weight import AffineWeight

        if not isinstance(weight, AffineWeight):
            # Fallback to Sage's action when possible.
            return self._w.action(weight)

        translated = self._group._translate_affine_weight(self._beta, weight)
        new_finite = self._w.action(translated.finite_part)
        return AffineWeight(weight.algebra, new_finite, translated.level, translated.grade)

    __call__ = action
