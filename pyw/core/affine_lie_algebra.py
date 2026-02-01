"""
Affine Lie Algebra Utilities (Di Francesco Convention)

This module provides core functions for affine Lie algebra computations,
following the notation and conventions from:

    Di Francesco, Mathieu, Sénéchal - "Conformal Field Theory"
    Chapter 14: Affine Lie Algebras

Key conventions from Di Francesco:
    - Affine weights: ̂λ = (λ; k_λ; n_λ) = (finite part; level; L0 eigenvalue)
    - Affine roots: ̂α = (α; 0; n) = α + nδ
    - Imaginary root: δ = (0; 0; 1), with (δ, δ) = 0
    - Extra simple root: α₀ = (-θ; 0; 1) = -θ + δ
    - Level: k = Σ a_i^∨ λ_i = (̂λ, δ) = ̂λ(̂k)
    - Scalar product: (̂λ, ̂μ) = (λ, μ) + k_λ n_μ + k_μ n_λ
    - Affine Weyl reflection: s_̂α ̂λ = ̂λ - (̂λ, ̂α^∨) ̂α
    - Translation: t_α∨(λ; k; n) = (λ + kα∨; k; n + [|λ|² - |λ + kα∨|²]/(2k))

Reference: Di Francesco et al., Chapter 14, Eqs. (14.23), (14.63), (14.64), (14.69)
"""

from typing import Union, Tuple, Dict, Any, Optional
from sage.all import RootSystem, CartanType, WeylGroup, QQ, ZZ, Integer, Rational


class AffineLieAlgebra:
    """
    Core utilities for affine Lie algebra computations following Di Francesco.

    This class provides convenient access to:
    - Marks (a_i) and comarks (a_i^∨)
    - Scalar products for finite and affine weights
    - Weyl reflections (finite and affine)
    - Translation operators t_α∨
    - Special elements: δ, θ, α₀, ρ̂

    Parameters
    ----------
    cartan_type : tuple or list
        Cartan type specification, e.g., ['A', 2, 1] for affine A₂

    Examples
    --------
    >>> from pyw.core.affine_lie_algebra import AffineLieAlgebra
    >>> ala = AffineLieAlgebra(['A', 2, 1])
    >>> marks = ala.get_marks()
    >>> comarks = ala.get_comarks()
    >>> print(f"Marks: {marks}, Comarks: {comarks}")
    """

    def __init__(self, cartan_type: Union[tuple, list]) -> None:
        """
        Initialize the AffineLieAlgebra utilities.

        Parameters
        ----------
        cartan_type : tuple or list
            Cartan type, e.g., ('A', 2) for finite A₂
            or ['A', 2, 1] for affine A₂^(1)
        """
        self._cartan_type = cartan_type

        # Initialize SageMath structures
        self._root_system = None
        self._cartan_type_obj = None
        self._weight_lattice = None
        self._root_lattice = None
        self._coroot_lattice = None
        self._ambient_space = None
        self._weyl_group = None

        # Cache for marks and comarks
        self._marks = None
        self._comarks = None

        # Cache for special elements
        self._delta = None
        self._theta = None
        self._rho_hat = None

        # Finite type information (for affine types)
        self._finite_type = None
        self._finite_root_system = None

        self._setup()

    def _setup(self) -> None:
        """Initialize SageMath structures."""
        self._root_system = RootSystem(self._cartan_type)
        self._cartan_type_obj = self._root_system.cartan_type()
        self._weight_lattice = self._root_system.weight_lattice()
        self._root_lattice = self._root_system.root_lattice()
        self._coroot_lattice = self._root_system.coroot_lattice()

        # Get ambient space for scalar products
        self._ambient_space = self._root_system.ambient_space()

        # For affine types, also get the finite type information
        if self.is_affine:
            finite_type_str = str(self._cartan_type_obj).split()[1]
            # Parse the finite type (e.g., ['A', 2] from "['A', 2, 1]")
            if isinstance(self._cartan_type, list) and len(self._cartan_type) == 3:
                self._finite_type = (self._cartan_type[0], self._cartan_type[1])
                self._finite_root_system = RootSystem(self._finite_type)

    # ==========================================================================
    # Type Detection
    # ==========================================================================

    @property
    def is_affine(self) -> bool:
        """Check if this is an affine Cartan type."""
        return self._cartan_type_obj.is_affine()

    @property
    def is_finite(self) -> bool:
        """Check if this is a finite Cartan type."""
        return not self.is_affine

    @property
    def rank(self) -> int:
        """Return the rank of the (finite) Lie algebra."""
        return self._cartan_type_obj.rank()

    @property
    def affine_rank(self) -> int:
        """Return the affine rank (rank + 1 for affine types)."""
        if self.is_affine:
            return self._cartan_type_obj.rank()
        return self.rank

    # ==========================================================================
    # Marks and Comarks (Di Francesco §14.1.3)
    # ==========================================================================

    def get_marks(self) -> Dict[int, int]:
        """
        Get the marks (a_i) for this Cartan type.

        Marks are the coefficients in the expansion of the highest root θ
        in terms of simple roots: θ = Σ a_i α_i.

        For affine types, a_0 = 1 by definition.

        Returns
        -------
        dict
            Dictionary mapping simple root indices to marks a_i

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> ala.get_marks()
        {0: 1, 1: 1, 2: 1}

        >>> ala = AffineLieAlgebra(['G', 2, 1])
        >>> ala.get_marks()
        {0: 1, 1: 3, 2: 1}

        Notes
        -----
        See Di Francesco Eq. (14.38) and Fig. 14.1 for marks values.
        """
        if self._marks is None:
            ct = self._cartan_type_obj
            # SageMath uses .c() for marks (a_i)
            sage_marks = ct.c()
            self._marks = {i: int(sage_marks[i]) for i in sage_marks.keys()}
        return self._marks

    def get_comarks(self) -> Dict[int, int]:
        """
        Get the comarks (a_i^∨) for this Cartan type.

        Comarks are the coefficients in the expansion of the highest
        coroot θ^∨ in terms of simple coroots.

        For affine types, a_0^∨ = 1 by definition.

        Returns
        -------
        dict
            Dictionary mapping simple root indices to comarks a_i^∨

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> ala.get_comarks()
        {0: 1, 1: 1, 2: 1}

        >>> ala = AffineLieAlgebra(['C', 2, 1])
        >>> ala.get_comarks()
        {0: 1, 1: 1, 2: 1}

        Notes
        -----
        See Di Francesco Eq. (14.39) and Fig. 14.1 for comarks values.
        For simply-laced algebras, marks = comarks.
        """
        if self._comarks is None:
            ct = self._cartan_type_obj
            # Comarks are marks of the dual Cartan type
            sage_comarks = ct.dual().c()
            self._comarks = {i: int(sage_comarks[i]) for i in sage_comarks.keys()}
        return self._comarks

    def dual_coxeter_number(self) -> int:
        """
        Get the dual Coxeter number g^∨ = Σ a_i^∨.

        Returns
        -------
        int
            The dual Coxeter number

        Examples
        --------
        >>> AffineLieAlgebra(['A', 2, 1]).dual_coxeter_number()
        3
        >>> AffineLieAlgebra(['B', 2, 1]).dual_coxeter_number()
        3
        >>> AffineLieAlgebra(['G', 2, 1]).dual_coxeter_number()
        4

        Notes
        -----
        See Di Francesco Eq. (14.42): g = Σ a_i^∨
        """
        return sum(self.get_comarks().values())

    # ==========================================================================
    # Scalar Products (Di Francesco §14.1.2)
    # ==========================================================================

    def scalar_product(self, lambda1, lambda2) -> Any:
        """
        Compute the scalar product (λ, μ) for two weights.

        This method handles both finite weights and AffineWeight objects:
        - For finite weights: computes (λ, μ) using Killing form
        - For AffineWeight: computes (̂λ, ̂μ) using Di Francesco Eq. (14.23)

        The Killing form normalization uses:
            (λ, μ) = Σ_{i,j} λ_i μ_j (C^{-1})_{ij}

        where λ_i, μ_j are the Dynkin labels and C is the Cartan matrix.

        Parameters
        ----------
        lambda1, lambda2
            Weights in the weight lattice, root lattice, or AffineWeight objects

        Returns
        -------
        The scalar product as a rational number

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2])
        >>> Lambda = ala.fundamental_weights()
        >>> ala.scalar_product(Lambda[1], Lambda[1])
        2/3
        >>> ala.scalar_product(Lambda[1], Lambda[2])
        1/3

        For affine weights (Di Francesco Eq. 14.23):

        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> Lambda_hat = ala.affine_fundamental_weights()
        >>> ala.scalar_product(Lambda_hat[1], Lambda_hat[2])
        1/3
        >>> ala.scalar_product(Lambda_hat[0], Lambda_hat[0])
        0

        Notes
        -----
        For AffineWeight objects, uses Di Francesco Eq. (14.23):
            (̂λ, ̂μ) = (λ, μ) + k_λ n_μ + k_μ n_λ

        The Killing form normalization differs from SageMath's ambient space
        Euclidean inner product. This method implements the correct Kac-Moody
        normalization where (Λ_i, Λ_j) = (C^{-1})_{ij}.
        """
        from sage.all import matrix, QQ

        # Check if inputs are AffineWeight objects
        from .affine_weight import AffineWeight

        if isinstance(lambda1, AffineWeight) and isinstance(lambda2, AffineWeight):
            # Use Di Francesco Eq. (14.23): (̂λ, ̂μ) = (λ, μ) + k_λ n_μ + k_μ n_λ
            finite_part = self._finite_scalar_product(lambda1.finite_part, lambda2.finite_part)
            cross_terms = lambda1.level * lambda2.grade + lambda2.level * lambda1.grade
            return finite_part + cross_terms
        elif isinstance(lambda1, AffineWeight) or isinstance(lambda2, AffineWeight):
            raise TypeError(
                "Both arguments must be AffineWeight objects, or both must be finite weights. "
                "Cannot mix AffineWeight with finite weights."
            )

        # Finite weights case
        return self._finite_scalar_product(lambda1, lambda2)

    def _finite_scalar_product(self, lambda1, lambda2) -> Any:
        """
        Compute the scalar product (λ, μ) for two finite weights.

        Internal method used by scalar_product().
        """
        from sage.all import matrix, QQ

        # Get the finite root system (for affine types)
        if self.is_affine:
            rs = self._finite_root_system
        else:
            rs = self._root_system

        # Get Cartan matrix and its inverse
        C = matrix(QQ, rs.cartan_matrix())
        C_inv = C.inverse()

        # Get Dynkin labels (coefficients of fundamental weights)
        # For a weight λ = Σ λ_i Λ_i, we have (λ, μ) = Σ_{i,j} λ_i μ_j (C^{-1})_{ij}

        # Extract Dynkin labels from monomial_coefficients
        # This gives us the coefficients in the fundamental weight basis
        mc1 = lambda1.monomial_coefficients()
        mc2 = lambda2.monomial_coefficients()

        # For affine types, we need to extract only the finite part (indices > 0)
        # For finite types, we use all indices
        if self.is_affine:
            # Get finite index set (excluding 0)
            finite_indices = [i for i in rs.index_set()]
            # Map affine indices to finite indices
            # For affine A2^(1): indices are 0,1,2 -> finite indices are 1,2
            # We need to map: affine index i -> finite index (i-1) for i > 0
            index_map = {idx: i for i, idx in enumerate(finite_indices)}
        else:
            finite_indices = list(rs.index_set())
            index_map = {idx: i for i, idx in enumerate(finite_indices)}

        # Build Dynkin label vectors
        n = C_inv.nrows()
        lambda1_labels = [0] * n
        lambda2_labels = [0] * n

        for idx, coeff in mc1.items():
            if self.is_affine and idx == 0:
                # Skip the affine index 0
                continue
            if idx in index_map:
                lambda1_labels[index_map[idx]] = coeff

        for idx, coeff in mc2.items():
            if self.is_affine and idx == 0:
                # Skip the affine index 0
                continue
            if idx in index_map:
                lambda2_labels[index_map[idx]] = coeff

        # Compute (λ, μ) = Σ_{i,j} λ_i μ_j (C^{-1})_{ij}
        result = 0
        for i in range(n):
            for j in range(n):
                result += lambda1_labels[i] * lambda2_labels[j] * C_inv[i, j]

        return result

    def affine_scalar_product(
        self, lambda_finite, k_lambda, n_lambda, mu_finite, k_mu, n_mu
    ) -> Any:
        """
        Compute the affine scalar product (̂λ, ̂μ).

        Following Di Francesco Eq. (14.23):
            (̂λ, ̂μ) = (λ, μ) + k_λ n_μ + k_μ n_λ

        where ̂λ = (λ; k_λ; n_λ), ̂μ = (μ; k_μ; n_μ)

        Parameters
        ----------
        lambda_finite : weight
            Finite part of ̂λ
        k_lambda : number
            Level of ̂λ (k_λ = ̂λ(̂k))
        n_lambda : number
            L₀ eigenvalue of ̂λ (with sign convention -L₀)
        mu_finite : weight
            Finite part of ̂μ
        k_mu : number
            Level of ̂μ (k_μ = ̂μ(̂k))
        n_mu : number
            L₀ eigenvalue of ̂μ

        Returns
        -------
        The affine scalar product (̂λ, ̂μ)

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> Lambda = ala.fundamental_weights()
        >>> # ̂λ = (Λ₁; 1; 0), ̂μ = (Λ₂; 1; 0)
        >>> ala.affine_scalar_product(Lambda[1], 1, 0, Lambda[2], 1, 0)
        -1/2

        Notes
        -----
        Di Francesco Eq. (14.23): (̂λ, ̂μ) = (λ, μ) + k_λ n_μ + k_μ n_λ

        Note: Di Francesco defines n_λ as ̂λ(-L₀), so our convention
        matches his: ̂λ = (λ; k_λ; n_λ) means ̂λ(-L₀) = n_λ.
        """
        finite_part = self._finite_scalar_product(lambda_finite, mu_finite)
        cross_terms = k_lambda * n_mu + k_mu * n_lambda
        return finite_part + cross_terms

    # ==========================================================================
    # Weyl Reflections (Di Francesco §14.1.6)
    # ==========================================================================

    def weyl_reflection(self, alpha, weight, alpha_vee=None) -> Any:
        """
        Compute the finite Weyl reflection s_α(λ).

        Following Di Francesco Eq. (14.63) (finite version):
            s_α(λ) = λ - (λ, α^∨) α

        Parameters
        ----------
        alpha : root
            The simple root α defining the reflection
        weight : weight
            The weight λ to reflect
        alpha_vee : coroot, optional
            The coroot α^∨. If None, uses weight.dot(alpha) or default.

        Returns
        -------
        The reflected weight s_α(λ)

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2])
        >>> Lambda = ala.fundamental_weights()
        >>> alpha = ala.simple_roots()
        >>> ala.weyl_reflection(alpha[1], Lambda[1])
        -Lambda[1] + Lambda[2]

        Notes
        -----
        Di Francesco Eq. (14.63) for finite case.
        The parameter order (alpha, weight) matches the mathematical notation s_α(λ).
        """
        # Use SageMath's built-in reflection if available
        try:
            # Try using weight lattice simple_reflection
            return weight.simple_reflection(self._root_index(alpha))
        except (AttributeError, TypeError):
            pass

        # Fallback: manual computation
        if alpha_vee is None:
            # Try to get coroot from root lattice
            try:
                alpha_vee = self._coroot_lattice.simple_coroot(self._root_index(alpha))
            except:
                # Use scalar product with the root itself
                pairing = self.scalar_product(weight, alpha)
                norm_sq = self.scalar_product(alpha, alpha)
                return weight - (2 * pairing / norm_sq) * alpha

        # Compute (λ, α^∨) and return λ - (λ, α^∨)α
        try:
            pairing = weight.dot(alpha_vee)
        except AttributeError:
            # Fallback to scalar product
            pairing = self.scalar_product(weight, alpha_vee)

        return weight - pairing * alpha

    def affine_weyl_reflection(
        self, alpha_hat: "AffineWeight", lambda_hat: "AffineWeight"
    ) -> "AffineWeight":
        """
        Compute the affine Weyl reflection s_̂α(̂λ).

        Following Di Francesco Eq. (14.64):
            s_̂α ̂λ = (λ - [(λ, α) + k m] α^∨; k;
                        n - [(λ, α) + k m] 2m/|α|²)

        where ̂λ = (λ; k; n), ̂α = (α; 0; m)

        Parameters
        ----------
        alpha_hat : AffineWeight
            The affine root ̂α = (α; 0; m) defining the reflection.
            Must have level = 0 (roots have zero level).
        lambda_hat : AffineWeight
            The affine weight ̂λ = (λ; k; n) to reflect.

        Returns
        -------
        AffineWeight
            The reflected affine weight s_̂α(̂λ)

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> alpha_hat = ala.affine_simple_roots()
        >>> Lambda_hat = ala.affine_fundamental_weights()
        >>> # Reflect ̂Λ₁ by ̂α₁
        >>> result = ala.affine_weyl_reflection(alpha_hat[1], Lambda_hat[1])
        >>> print(result)

        Notes
        -----
        Di Francesco Eq. (14.64):
            s_̂α ̂λ = (λ - [(λ, α) + k m] α^∨;
                       k;
                       n - [(λ, α) + k m] 2m/|α|²)

        For m = 0, this reduces to the finite Weyl reflection.
        The parameter order (alpha_hat, lambda_hat) matches s_̂α(̂λ).
        """
        from .affine_weight import AffineWeight

        # Extract components from AffineWeight objects
        alpha_finite = alpha_hat.finite_part
        m_alpha = alpha_hat.grade  # grade is the m component for roots

        lambda_finite = lambda_hat.finite_part
        k_lambda = lambda_hat.level
        n_lambda = lambda_hat.grade

        # Validate that alpha_hat is a root (level = 0)
        if alpha_hat.level != 0:
            raise ValueError(
                f"alpha_hat must be an affine root with level=0, got level={alpha_hat.level}"
            )

        # Compute (λ, α) - finite scalar product
        lambda_dot_alpha = self._finite_scalar_product(lambda_finite, alpha_finite)

        # Get |α|² for normalization
        alpha_norm_sq = self._finite_scalar_product(alpha_finite, alpha_finite)

        # Compute (λ, α^∨) = 2(λ, α)/|α|²
        if alpha_norm_sq == 0:
            lambda_dot_alpha_vee = 0
        else:
            lambda_dot_alpha_vee = 2 * lambda_dot_alpha / alpha_norm_sq

        # Full coefficient for Di Francesco Eq. (14.64):
        # [(λ, α) + k*m] for the finite part
        full_coeff = lambda_dot_alpha + k_lambda * m_alpha

        # New finite part: s_α(λ) for m=0, or more complex for m≠0
        # Use SageMath's simple_reflection when possible
        try:
            # Try to use weight's simple_reflection method
            # This works when alpha is a simple root and lambda is in weight space
            new_lambda = lambda_finite.simple_reflection(self._finite_root_index(alpha_finite))
        except (AttributeError, TypeError, ValueError):
            # Fallback: compute manually
            # s_α(λ) = λ - (λ, α^∨) α
            # Need to express α in the same space as λ
            # For weight space, α_i = Σ_j C_ij Λ_j
            new_lambda = self._apply_reflection_manually(
                lambda_finite, alpha_finite, lambda_dot_alpha_vee
            )

        # Level is unchanged
        new_k = k_lambda

        # New L₀ eigenvalue: n - [(λ, α) + k*m] * 2m/|α|²
        if alpha_norm_sq == 0 or m_alpha == 0:
            n_correction = 0
        else:
            n_correction = full_coeff * (2 * m_alpha / alpha_norm_sq)
        new_n = n_lambda - n_correction

        return AffineWeight(self, new_lambda, level=new_k, grade=new_n)

    def _finite_root_index(self, root) -> int:
        """Get the index of a simple root in the finite root system."""
        if self.is_affine:
            alpha = self._finite_root_system.root_lattice().simple_roots()
        else:
            alpha = self._root_lattice.simple_roots()

        for i, a in alpha.items():
            if root == a:
                return i
        raise ValueError("Root not found in finite simple roots")

    def _apply_reflection_manually(self, weight, alpha, pairing):
        """Apply Weyl reflection manually when simple_reflection is unavailable."""
        from sage.all import matrix, QQ

        # Get the appropriate root system
        if self.is_affine:
            rs = self._finite_root_system
        else:
            rs = self._root_system

        ws = rs.weight_space()
        C = matrix(QQ, rs.cartan_matrix())
        finite_indices = list(rs.index_set())

        # Get alpha's coefficients in simple root basis
        alpha_coeffs = alpha.monomial_coefficients()

        # Build alpha in weight space: α_i = Σ_j C_ij Λ_j
        Lambda_ws = ws.fundamental_weights()
        alpha_in_ws = ws.zero()
        for root_idx, root_coeff in alpha_coeffs.items():
            i_pos = finite_indices.index(root_idx)
            for j_pos, j_idx in enumerate(finite_indices):
                alpha_in_ws += root_coeff * C[i_pos, j_pos] * Lambda_ws[j_idx]

        return weight - pairing * alpha_in_ws

    def simple_reflection(self, i: int, weight, k=None, n=None):
        """
        Apply simple reflection s_i to a weight.

        For affine types with k and n provided, computes affine reflection.
        For finite types or when k,n are None, computes finite reflection.

        Parameters
        ----------
        i : int
            Index of the simple root (0 for affine node in affine types)
        weight : weight
            The weight to reflect
        k : int, optional
            Level (for affine reflections)
        n : int, optional
            L₀ eigenvalue (for affine reflections)

        Returns
        -------
        The reflected weight, or tuple (weight, k, n) for affine case

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> Lambda = ala.fundamental_weights()
        >>> # Finite reflection s₁ on Λ₁
        >>> ala.simple_reflection(1, Lambda[1])
        -Lambda[1] + Lambda[2]
        """
        alpha = self.simple_roots().get(i)

        if k is None or n is None:
            # Finite reflection
            return self.weyl_reflection(alpha, weight)
        else:
            # Affine reflection with m = 0 (simple roots have m=0 for i > 0)
            # For i = 0, the affine simple root has m = 1
            from .affine_weight import AffineWeight

            # Get the affine simple root
            alpha_hat = AffineWeight.affine_simple_root(self, i)
            # Create affine weight from components
            lambda_hat = AffineWeight(self, weight, level=k, grade=n)
            # Perform affine reflection
            result = self.affine_weyl_reflection(alpha_hat, lambda_hat)
            # Return as tuple for backward compatibility
            return (result.finite_part, result.level, result.grade)

    # ==========================================================================
    # Translation Operators (Di Francesco §14.1.6)
    # ==========================================================================

    def translation(self, alpha_vee, weight: "AffineWeight") -> "AffineWeight":
        """
        Compute the translation t_α∨(̂λ) following Di Francesco Eq. (14.69).

        The translation operator t_α∨ acts on affine weights as:
            t_α∨(λ; k; n) = (λ + k α∨; k; n - [(λ, α∨) + k|α∨|²/2])

        This is the simplified form of the original formula:
            t_α∨(λ; k; n) = (λ + k α∨; k; n + [|λ|² - |λ + kα∨|²]/(2k))

        The simplified form is well-defined even when k = 0.

        Parameters
        ----------
        alpha_vee : coroot
            The coroot α^∨ defining the translation direction
        weight : AffineWeight
            The affine weight ̂λ = (λ; k; n)

        Returns
        -------
        AffineWeight
            The translated affine weight t_α∨(̂λ)

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> alpha_vee = ala.simple_coroots()[1]
        >>> w = AffineWeight(ala, ala.weight_lattice().zero(), level=1, grade=0)
        >>> result = ala.translation(alpha_vee, w)
        >>> # Result should be (α₁^∨; 1; -|α₁^∨|²/2)

        Notes
        -----
        Di Francesco Eq. (14.69) can be simplified:
            n' = n + [|λ|² - |λ + kα∨|²]/(2k)
               = n - (λ, α∨) - k|α∨|²/2

        This form avoids division by k and is valid for all k including k = 0.

        This is the Di Francesco convention. Note that Kac-Wakimoto
        uses the opposite sign: t_β^KW = t_{-β}^DF.

        The parameter order (coroot, weight) matches the mathematical
        notation t_α∨(λ).
        """
        from .affine_weight import AffineWeight

        lambda_finite = weight.finite_part
        k = weight.level
        n = weight.grade

        # Convert alpha_vee to the same space as lambda_finite if needed
        # Coroots are in coroot_lattice, but we need to add to weight_lattice/space
        alpha_vee_coeffs = alpha_vee.monomial_coefficients()
        if hasattr(lambda_finite, "parent"):
            parent = lambda_finite.parent()
            # Reconstruct alpha_vee in the weight space using simple roots
            simple_elements = (
                parent.simple_roots()
                if hasattr(parent, "simple_roots")
                else parent.fundamental_weights()
            )
            alpha_vee_in_parent = parent.zero()
            for idx, coeff in alpha_vee_coeffs.items():
                if idx in simple_elements.keys():
                    alpha_vee_in_parent += coeff * simple_elements[idx]
        else:
            alpha_vee_in_parent = alpha_vee

        # New finite part: λ + k α∨
        new_lambda = lambda_finite + k * alpha_vee_in_parent

        # Simplified grade correction (valid for all k including k = 0):
        # n' = n - (λ, α∨) - k|α∨|²/2
        lambda_dot_alpha_vee = self.scalar_product(lambda_finite, alpha_vee)
        alpha_vee_norm_sq = self.scalar_product(alpha_vee, alpha_vee)
        new_n = n - lambda_dot_alpha_vee - k * alpha_vee_norm_sq / 2

        return AffineWeight(self, new_lambda, k, new_n)

    def translation_by_root(self, alpha, weight: "AffineWeight") -> "AffineWeight":
        """
        Compute translation t_α where α is a root (not coroot).

        This converts the root to coroot using α^∨ = (2/|α|²) α.

        Parameters
        ----------
        alpha : root
            The root α defining the translation direction
        weight : AffineWeight
            The affine weight ̂λ = (λ; k; n)

        Returns
        -------
        AffineWeight
            The translated affine weight t_α(̂λ)

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> alpha = ala.simple_roots()[1]
        >>> w = AffineWeight(ala, ala.weight_space().zero(), level=1, grade=0)
        >>> result = ala.translation_by_root(alpha, w)

        Notes
        -----
        The parameter order (root, weight) matches the mathematical
        notation t_α(λ).
        """
        # Convert root to coroot: α^∨ = (2/|α|²) α
        alpha_norm_sq = self.scalar_product(alpha, alpha)
        if alpha_norm_sq == 0:
            raise ValueError("Cannot translate by imaginary root (zero norm)")

        alpha_vee = (2 / alpha_norm_sq) * alpha
        return self.translation(alpha_vee, weight)

    # ==========================================================================
    # Special Elements (δ, θ, α₀, ρ̂)
    # ==========================================================================

    def delta(self):
        """
        Get the imaginary root δ = (0; 0; 1).

        In Di Francesco's notation, δ has zero norm: (δ, δ) = 0.
        It generates the imaginary roots nδ.

        Returns
        -------
        A representation of δ (as a tuple or symbolic element)

        Notes
        -----
        Di Francesco Eq. (14.27): δ = (0; 0; 1)
        Di Francesco Eq. (14.31): (δ, δ) = 0
        """
        return (0, 0, 1)  # (finite; k-grade; n-grade)

    def theta(self):
        """
        Get the highest root θ of the finite Lie algebra.

        Returns
        -------
        The highest root θ

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> theta = ala.theta()
        >>> # For A₂, θ = α₁ + α₂

        Notes
        -----
        The highest root θ is used to define α₀ = -θ + δ.
        See Di Francesco Eq. (14.32).
        """
        if self._finite_root_system is not None:
            # Use finite root system
            Q_finite = self._finite_root_system.root_lattice()
            return Q_finite.highest_root()
        else:
            # Finite type - use current root system
            return self._root_lattice.highest_root()

    def alpha_0(self):
        """
        Get the extra simple root α₀ = (-θ; 0; 1) = -θ + δ.

        Returns
        -------
        The affine simple root α₀

        Notes
        -----
        Di Francesco Eq. (14.32): α₀ = (-θ; 0; 1) = -θ + δ
        """
        theta = self.theta()
        # α₀ = -θ + δ, but as a finite root it's -θ
        return -theta

    def rho_hat(self):
        """
        Get the affine Weyl vector ρ̂.

        Returns
        -------
        The affine Weyl vector ρ̂

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> rho_hat = ala.rho_hat()
        >>> # In Dynkin labels: ρ̂ = [1, 1, 1]

        Notes
        -----
        Di Francesco Eq. (14.62): ρ̂ = Σ ω̂_i = [1, 1, ..., 1]
        Note: ρ̂(k) = g (the dual Coxeter number)
        """
        if self._rho_hat is None:
            # For affine types, ρ̂ is sum of fundamental weights
            Lambda = self.fundamental_weights()
            self._rho_hat = sum(Lambda.values())
        return self._rho_hat

    # ==========================================================================
    # Access to SageMath Structures
    # ==========================================================================

    def root_system(self):
        """Get the underlying SageMath RootSystem."""
        return self._root_system

    def weight_lattice(self):
        """Get the weight lattice."""
        return self._weight_lattice

    def root_lattice(self):
        """Get the root lattice."""
        return self._root_lattice

    def coroot_lattice(self):
        """Get the coroot lattice."""
        return self._coroot_lattice

    def ambient_space(self):
        """Get the ambient space (for scalar products)."""
        return self._ambient_space

    def fundamental_weights(self, lattice: bool = False):
        """
        Get the fundamental weights Λ_i.

        Parameters
        ----------
        lattice : bool, optional
            If True, return weights from the weight lattice (integer coefficients only).
            If False (default), return weights from the weight space (rational coefficients).

        Returns
        -------
        dict mapping indices to fundamental weights

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> Lambda = ala.fundamental_weights()
        >>> Lambda[0]  # Basic fundamental weight ω̂₀ = (0; 1; 0)
        >>> Lambda[1] / 3  # Rational coefficients supported
        """
        if lattice:
            Lambda = self._weight_lattice.fundamental_weights()
        else:
            # Use weight space for rational coefficient support
            ws = self._root_system.weight_space()
            Lambda = ws.fundamental_weights()
        return {i: Lambda[i] for i in self._weight_lattice.index_set()}

    def simple_roots(self):
        """
        Get the simple roots α_i.

        Returns
        -------
        dict mapping indices to simple roots

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> alpha = ala.simple_roots()
        >>> alpha[0]  # α₀ = -θ + δ
        """
        alpha = self._root_lattice.simple_roots()
        return {i: alpha[i] for i in self._root_lattice.index_set()}

    def simple_coroots(self):
        """
        Get the simple coroots α_i^∨.

        Returns
        -------
        dict mapping indices to simple coroots

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> alpha_vee = ala.simple_coroots()
        >>> alpha_vee[0]  # α₀^∨
        """
        alpha_vee = self._coroot_lattice.simple_roots()
        return {i: alpha_vee[i] for i in self._coroot_lattice.index_set()}

    def cartan_matrix(self):
        """Get the Cartan matrix."""
        return self._cartan_type_obj.cartan_matrix()

    # ==========================================================================
    # Helper Methods
    # ==========================================================================

    def _root_index(self, root) -> int:
        """Get the index of a simple root."""
        alpha = self.simple_roots()
        for i, a in alpha.items():
            if root == a:
                return i
        raise ValueError("Root not found in simple roots")

    def level_from_dynkin(self, dynkin_labels: Dict[int, int]) -> int:
        """
        Compute the level k from Dynkin labels [λ₀, λ₁, ..., λ_r].

        Following Di Francesco Eq. (14.54):
            k = Σ a_i^∨ λ_i

        Parameters
        ----------
        dynkin_labels : dict
            Dictionary mapping indices i to Dynkin labels λ_i

        Returns
        -------
        int
            The level k

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> # Weight [1, 0, 0] has level k = 1*1 + 1*0 + 1*0 = 1
        >>> ala.level_from_dynkin({0: 1, 1: 0, 2: 0})
        1

        Notes
        -----
        Di Francesco Eq. (14.54): k = Σ a_i^∨ λ_i
        """
        comarks = self.get_comarks()
        k = sum(comarks[i] * dynkin_labels.get(i, 0) for i in comarks.keys())
        return k

    def lambda0_from_level(self, k: int, finite_labels: Dict[int, int]) -> int:
        """
        Compute λ₀ from level k and finite Dynkin labels.

        Following Di Francesco Eq. (14.57):
            λ₀ = k - (λ, θ) = k - Σ a_i λ_i

        Parameters
        ----------
        k : int
            The level
        finite_labels : dict
            Dictionary of finite Dynkin labels (indices 1, 2, ..., r)

        Returns
        -------
        int
            The zeroth Dynkin label λ₀

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> # For k=2 and finite labels {1: 1, 2: 0}:
        >>> # λ₀ = 2 - (1*1 + 1*0) = 1
        >>> ala.lambda0_from_level(2, {1: 1, 2: 0})
        1

        Notes
        -----
        Di Francesco Eq. (14.57): λ₀ = k - (λ, θ)
        """
        marks = self.get_marks()
        # Sum over finite marks: (λ, θ) = Σ_{i>0} a_i λ_i
        theta_dot_lambda = sum(marks[i] * finite_labels.get(i, 0) for i in marks.keys() if i != 0)
        return k - theta_dot_lambda

    def is_dominant(self, dynkin_labels: Dict[int, int]) -> bool:
        """
        Check if a weight is dominant (all Dynkin labels ≥ 0).

        Parameters
        ----------
        dynkin_labels : dict
            Dictionary of Dynkin labels

        Returns
        -------
        bool
            True if dominant, False otherwise

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> ala.is_dominant({0: 1, 1: 0, 2: 0})
        True
        >>> ala.is_dominant({0: -1, 1: 1, 2: 1})
        False
        """
        return all(l >= 0 for l in dynkin_labels.values())

    # ==========================================================================
    # Di Francesco Notation Factory Methods
    # ==========================================================================

    def affine_weight(self, finite_part, level=0, grade=0):
        """
        Create an AffineWeight in Di Francesco notation: (λ; k; n).

        Parameters
        ----------
        finite_part : weight
            The finite part λ
        level : number, optional
            The level k (default: 0)
        grade : number, optional
            The L₀ eigenvalue n (default: 0)

        Returns
        -------
        AffineWeight
            The affine weight (λ; k; n)

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> ws = ala._finite_root_system.weight_space()
        >>> Lambda = ws.fundamental_weights()
        >>> w = ala.affine_weight(Lambda[1], level=1, grade=0)
        """
        from .affine_weight import AffineWeight

        return AffineWeight(self, finite_part, level, grade)

    def affine_fundamental_weights(self):
        """
        Get affine fundamental weights in Di Francesco notation.

        Returns ̂Λ_i = (Λ_i; a_i^∨; 0) for i = 0, 1, ..., r.

        For i = 0: ̂Λ₀ = (0; 1; 0)
        For i > 0: ̂Λ_i = (Λ_i; 1; 0)

        Returns
        -------
        dict
            Dictionary mapping indices to AffineWeight objects

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> Lambda_hat = ala.affine_fundamental_weights()
        >>> print(Lambda_hat[0])  # (0; 1; 0)
        >>> print(Lambda_hat[1])  # (Lambda[1]; 1; 0)
        """
        from .affine_weight import AffineWeight

        if not self.is_affine:
            raise ValueError("Algebra must be affine type")

        index_set = list(self._weight_lattice.index_set())
        return {i: AffineWeight.affine_fundamental_weight(self, i) for i in index_set}

    def affine_simple_roots(self):
        """
        Get affine simple roots in Di Francesco notation.

        Returns ̂α_i for i = 0, 1, ..., r.

        For i > 0: ̂α_i = (α_i; 0; 0)
        For i = 0: ̂α₀ = (-θ; 0; 1) = -θ + δ

        Returns
        -------
        dict
            Dictionary mapping indices to AffineWeight objects

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> alpha_hat = ala.affine_simple_roots()
        >>> print(alpha_hat[0])  # (-θ; 0; 1)
        >>> print(alpha_hat[1])  # (alpha[1]; 0; 0)
        """
        from .affine_weight import AffineWeight

        if not self.is_affine:
            raise ValueError("Algebra must be affine type")

        index_set = list(self._root_lattice.index_set())
        return {i: AffineWeight.affine_simple_root(self, i) for i in index_set}

    def affine_delta(self):
        """
        Get the imaginary root δ = (0; 0; 1) in Di Francesco notation.

        The imaginary root satisfies (δ, δ) = 0 and (δ, ̂λ) = k_λ.

        Returns
        -------
        AffineWeight
            The imaginary root δ

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> delta = ala.affine_delta()
        >>> delta.norm_squared()  # 0
        """
        from .affine_weight import AffineWeight

        return AffineWeight.delta(self)

    def affine_rho(self):
        """
        Get the affine Weyl vector ρ̂ = (ρ; g; 0) in Di Francesco notation.

        The level of ρ̂ equals the dual Coxeter number g.

        Returns
        -------
        AffineWeight
            The affine Weyl vector ρ̂

        Examples
        --------
        >>> ala = AffineLieAlgebra(['A', 2, 1])
        >>> rho_hat = ala.affine_rho()
        >>> rho_hat.level  # 3 (dual Coxeter number for A₂)
        """
        from .affine_weight import AffineWeight

        return AffineWeight.rho_hat(self)

    def __repr__(self) -> str:
        return f"AffineLieAlgebra({self._cartan_type})"


# =============================================================================
# Convenience Functions (non-class interface)
# =============================================================================


def scalar_product(lambda1, lambda2, cartan_type=None):
    """
    Compute the scalar product (λ, μ) for two weights.

    This is a convenience function that creates an AffineLieAlgebra
    instance if needed.

    Parameters
    ----------
    lambda1, lambda2 : weights
        Weights from the same root system
    cartan_type : tuple or list, optional
        Cartan type if not inferrable from the weights

    Returns
    -------
    The scalar product (λ, μ)

    Examples
    --------
    >>> from pyw.core.affine_lie_algebra import scalar_product
    >>> from sage.all import RootSystem
    >>> R = RootSystem(['A', 2])
    >>> Lambda = R.weight_lattice().fundamental_weights()
    >>> scalar_product(Lambda[1], Lambda[2])
    -1/2
    """
    if cartan_type is not None:
        ala = AffineLieAlgebra(cartan_type)
    else:
        # Try to infer from weights
        # This is a simplified version - in practice you'd need
        # to extract the cartan type from the weight objects
        raise NotImplementedError("cartan_type must be specified")
    return ala.scalar_product(lambda1, lambda2)


def get_marks(cartan_type: Union[tuple, list]) -> Dict[int, int]:
    """
    Get the marks (a_i) for a Cartan type.

    Convenience function for AffineLieAlgebra.get_marks().

    Parameters
    ----------
    cartan_type : tuple or list
        Cartan type specification

    Returns
    -------
    dict
        Dictionary mapping indices to marks

    Examples
    --------
    >>> from pyw.core.affine_lie_algebra import get_marks
    >>> get_marks(['A', 2, 1])
    {0: 1, 1: 1, 2: 1}
    >>> get_marks(['G', 2, 1])
    {0: 1, 1: 3, 2: 1}
    """
    ala = AffineLieAlgebra(cartan_type)
    return ala.get_marks()


def get_comarks(cartan_type: Union[tuple, list]) -> Dict[int, int]:
    """
    Get the comarks (a_i^∨) for a Cartan type.

    Convenience function for AffineLieAlgebra.get_comarks().

    Parameters
    ----------
    cartan_type : tuple or list
        Cartan type specification

    Returns
    -------
    dict
        Dictionary mapping indices to comarks

    Examples
    --------
    >>> from pyw.core.affine_lie_algebra import get_comarks
    >>> get_comarks(['A', 2, 1])
    {0: 1, 1: 1, 2: 1}
    >>> get_comarks(['F', 4, 1])
    {0: 1, 1: 1, 2: 1, 3: 2, 4: 2}
    """
    ala = AffineLieAlgebra(cartan_type)
    return ala.get_comarks()


def weyl_reflection(alpha, weight, cartan_type: Union[tuple, list]):
    """
    Compute the Weyl reflection s_α(λ).

    Convenience function for AffineLieAlgebra.weyl_reflection().

    Parameters
    ----------
    alpha : root
        The root defining the reflection
    weight : weight
        The weight to reflect
    cartan_type : tuple or list
        Cartan type specification

    Returns
    -------
    The reflected weight s_α(λ)

    Examples
    --------
    >>> from pyw.core.affine_lie_algebra import weyl_reflection
    >>> from sage.all import RootSystem
    >>> R = RootSystem(['A', 2])
    >>> Lambda = R.weight_lattice().fundamental_weights()
    >>> alpha = R.root_lattice().simple_roots()
    >>> weyl_reflection(alpha[1], Lambda[1], ['A', 2])
    -Lambda[1] + Lambda[2]

    Notes
    -----
    The parameter order (alpha, weight) matches the mathematical notation s_α(λ).
    """
    ala = AffineLieAlgebra(cartan_type)
    return ala.weyl_reflection(alpha, weight)
