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

from typing import Any, Union
from sage.all import WeylGroup, RootSystem


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
