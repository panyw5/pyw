"""
Kazhdan-Lusztig Polynomials for Coxeter Groups

This module provides computation of Kazhdan-Lusztig polynomials and their
inverses, essential for character formulas of admissible representations.

Key features:
    - Standard KL polynomials P_{x,y}(q)
    - Inverse KL polynomials Q̃_{x,y}(q)
    - Parabolic (coset) KL polynomials
    - Integration with SageMath and coxeter3
    - Caching for expensive computations

The inverse KL polynomial Q̃_{x,y}(1) appears in the Kazhdan-Lusztig
character formula for admissible modules:

    ch(L_λ) = Σ Q̃_{[w],[w']}(1) · ch(M(w'·Λ))

References:
    - Kazhdan, D., Lusztig, G. "Representations of Coxeter groups..."
    - Soergel, W. "Kazhdan-Lusztig polynomials and a combinatoric for tilting modules"
    - Cordova, Gaiotto, Shao "Infrared Computations of Defect Schur Indices" (Eq. C.28)
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Union

from sage.all import QQ, SR, WeylGroup, matrix, var

if TYPE_CHECKING:
    from .bruhat import BruhatOrder, CosetRepresentative, ParabolicSubgroup


class KazhdanLusztigPolynomials:
    """
    Compute Kazhdan-Lusztig polynomials for Coxeter/Weyl groups.

    This class provides access to both standard KL polynomials P_{x,y}(q)
    and inverse KL polynomials Q̃_{x,y}(q), with support for parabolic
    (coset) versions needed for admissible character computations.

    The implementation uses multiple backends:
    1. SageMath's native KazhdanLusztigPolynomial (default)
    2. coxeter3 via invpol method (if available, faster for large groups)

    Parameters
    ----------
    coxeter_group : WeylGroup or CoxeterGroup
        The Coxeter group
    cache_dir : Path, optional
        Directory for caching computed polynomials

    Examples
    --------
    >>> from sage.all import WeylGroup
    >>> W = WeylGroup(['A', 2])
    >>> kl = KazhdanLusztigPolynomials(W)
    >>> s1 = W.simple_reflection(1)
    >>> s2 = W.simple_reflection(2)
    >>> kl.P(W.one(), s1 * s2)  # P_{e, s1s2}(q)
    1

    Notes
    -----
    The inverse KL polynomial Q̃_{x,y}(q) is related to P_{x,y}(q) by:
        Σ_z P_{x,z}(q) Q̃_{z,y}(q) = δ_{x,y}

    For the character formula, we need Q̃_{x,y}(1).
    """

    def __init__(
        self,
        coxeter_group: Any,
        cache_dir: Optional[Path] = None,
    ) -> None:
        """
        Initialize KL polynomial calculator.

        Parameters
        ----------
        coxeter_group : WeylGroup or CoxeterGroup
            The Coxeter group
        cache_dir : Path, optional
            Directory for caching (default: ~/.pyw/kl_cache)
        """
        self._W = coxeter_group
        self._cartan_type = coxeter_group.cartan_type()

        # Setup caching
        if cache_dir is None:
            cache_dir = Path.home() / ".pyw" / "kl_cache"
        self._cache_dir = cache_dir
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        # In-memory cache
        self._P_cache: Dict[Tuple[Any, Any], Any] = {}
        self._Q_cache: Dict[Tuple[Any, Any], Any] = {}
        self._invpol_cache: Dict[Tuple[Any, Any], Any] = {}

        # Initialize SageMath KL calculator
        self._sage_kl = None
        self._coxeter3 = None

        self._setup_backends()

    def _setup_backends(self) -> None:
        """Initialize computation backends."""
        # Try SageMath's native KL
        try:
            from sage.combinat.kazhdan_lusztig import KazhdanLusztigPolynomial
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

            # Use polynomial ring instead of symbolic variable to avoid subs() issues
            R = PolynomialRing(QQ, "q")
            q = R.gen()
            self._sage_kl = KazhdanLusztigPolynomial(self._W, q)
            self._q = q  # Store for substitution
        except ImportError:
            pass

        # Try coxeter3 backend (faster for large groups)
        try:
            from coxeter3_sage import Coxeter3

            # Convert Cartan type to coxeter3 format
            ct = self._cartan_type
            if hasattr(ct, "letter"):
                letter = ct.letter()
                rank = ct.rank()
                if ct.is_affine():
                    # Affine type: e.g., "A~2" for A_2^(1)
                    coxeter3_type = f"{letter}~{rank - 1}"
                else:
                    coxeter3_type = f"{letter}{rank}"
                self._coxeter3 = Coxeter3(coxeter3_type)
        except (ImportError, Exception):
            # coxeter3 not available or type not supported
            pass

    @property
    def weyl_group(self) -> Any:
        """The underlying Weyl/Coxeter group."""
        return self._W

    @property
    def cartan_type(self) -> Any:
        """The Cartan type."""
        return self._cartan_type

    # =========================================================================
    # Standard KL Polynomials P_{x,y}(q)
    # =========================================================================

    def P(self, x: Any, y: Any, at_one: bool = False) -> Any:
        """
        Compute the Kazhdan-Lusztig polynomial P_{x,y}(q).

        P_{x,y}(q) is defined for x ≤ y in Bruhat order and satisfies:
        - P_{x,x}(q) = 1
        - deg(P_{x,y}) ≤ (ℓ(y) - ℓ(x) - 1) / 2

        Parameters
        ----------
        x, y : Weyl group elements
            Elements with x ≤ y in Bruhat order
        at_one : bool
            If True, return P_{x,y}(1) instead of the polynomial

        Returns
        -------
        polynomial or int
            P_{x,y}(q) or P_{x,y}(1)

        Examples
        --------
        >>> W = WeylGroup(['A', 2])
        >>> kl = KazhdanLusztigPolynomials(W)
        >>> kl.P(W.one(), W.long_element())
        1
        """
        x = self._ensure_element(x)
        y = self._ensure_element(y)

        # Check cache
        cache_key = (self._element_key(x), self._element_key(y))
        if cache_key in self._P_cache:
            p = self._P_cache[cache_key]
            return p.subs({self._q: 1}) if at_one else p

        # Compute using SageMath
        if self._sage_kl is not None:
            p = self._sage_kl.P(x, y)
            self._P_cache[cache_key] = p
            return p.subs({self._q: 1}) if at_one else p

        # Fallback: P_{x,y} = 1 for x = y, need recursive formula otherwise
        raise NotImplementedError(
            "KL polynomial computation requires SageMath's KazhdanLusztigPolynomial"
        )

    # =========================================================================
    # Inverse KL Polynomials Q̃_{x,y}(q)
    # =========================================================================

    def Q_tilde(self, x: Any, y: Any, at_one: bool = True) -> Any:
        """
        Compute the inverse Kazhdan-Lusztig polynomial Q̃_{x,y}(q).

        Q̃_{x,y}(q) satisfies:
            Σ_z P_{x,z}(q) Q̃_{z,y}(q) = δ_{x,y}

        This is the polynomial needed for the KL character formula.

        Parameters
        ----------
        x, y : Weyl group elements
            Elements with x ≤ y in Bruhat order
        at_one : bool
            If True (default), return Q̃_{x,y}(1)

        Returns
        -------
        polynomial or int
            Q̃_{x,y}(q) or Q̃_{x,y}(1)

        Examples
        --------
        >>> W = WeylGroup(['A', 2])
        >>> kl = KazhdanLusztigPolynomials(W)
        >>> kl.Q_tilde(W.one(), W.long_element(), at_one=True)
        1
        """
        x = self._ensure_element(x)
        y = self._ensure_element(y)

        # Check cache
        cache_key = (self._element_key(x), self._element_key(y))
        if cache_key in self._Q_cache:
            q_val = self._Q_cache[cache_key]
            if at_one and hasattr(q_val, "subs"):
                return q_val.subs({self._q: 1})
            return q_val

        # Try coxeter3's invpol (fastest)
        if self._coxeter3 is not None and hasattr(self._coxeter3, "invpol"):
            try:
                result = self._coxeter3.invpol(x, y)
                self._Q_cache[cache_key] = result
                if at_one and hasattr(result, "subs"):
                    return result.subs({self._q: 1})
                return result
            except Exception:
                pass

        # Fallback: compute by inverting the P-matrix on the interval
        result = self._compute_Q_tilde_by_inversion(x, y)
        self._Q_cache[cache_key] = result

        if at_one and hasattr(result, "subs"):
            return result.subs({self._q: 1})
        return result

    def invpol(self, x: Any, y: Any) -> Any:
        """
        Compute inverse KL polynomial using coxeter3's invpol.

        This is a direct interface to coxeter3's invpol command,
        which computes Q̃_{x,y}(q) efficiently.

        Parameters
        ----------
        x, y : Weyl group elements
            Elements with x ≤ y in Bruhat order

        Returns
        -------
        polynomial
            Q̃_{x,y}(q) as a symbolic expression

        Raises
        ------
        RuntimeError
            If coxeter3 is not available
        """
        if self._coxeter3 is None:
            raise RuntimeError(
                "coxeter3 backend not available. "
                "Install coxeter3_sage and ensure invpol method is added."
            )

        x = self._ensure_element(x)
        y = self._ensure_element(y)

        # Check cache
        cache_key = (self._element_key(x), self._element_key(y))
        if cache_key in self._invpol_cache:
            return self._invpol_cache[cache_key]

        # Call coxeter3
        result = self._coxeter3.invpol(x, y)
        self._invpol_cache[cache_key] = result
        return result

    def _compute_Q_tilde_by_inversion(self, x: Any, y: Any) -> Any:
        """
        Compute Q̃_{x,y}(1) by inverting the P-matrix on [x, y].

        The inverse KL polynomials satisfy:
            Σ_z P_{x,z} Q̃_{z,y} = δ_{x,y}

        So Q̃ = P^{-1} on the Bruhat interval.
        """
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self._W)

        # Get the Bruhat interval [x, y]
        if not bruhat.le(x, y):
            return 0

        interval = bruhat.interval(x, y)
        n = len(interval)

        if n == 1:
            # x = y
            return 1

        # Build the P-matrix: P[i,j] = P_{interval[i], interval[j]}(1)
        # Only upper triangular (i ≤ j in Bruhat order)
        P_mat = matrix(QQ, n, n)
        for i, w_i in enumerate(interval):
            for j, w_j in enumerate(interval):
                if bruhat.le(w_i, w_j):
                    P_mat[i, j] = self.P(w_i, w_j, at_one=True)

        # Invert to get Q̃-matrix
        try:
            Q_mat = P_mat.inverse()
        except Exception as e:
            raise RuntimeError(f"Failed to invert P-matrix: {e}")

        # Find indices of x and y in the interval
        x_idx = interval.index(x)
        y_idx = interval.index(y)

        return Q_mat[x_idx, y_idx]

    # =========================================================================
    # Parabolic (Coset) KL Polynomials
    # =========================================================================

    def parabolic_Q_tilde(
        self,
        coset_x: "CosetRepresentative",
        coset_y: "CosetRepresentative",
        at_one: bool = True,
    ) -> Any:
        """
        Compute parabolic inverse KL polynomial Q̃_{[x],[y]}(q).

        For cosets [x], [y] in W / W_I (or W_I \\ W), the parabolic
        inverse KL polynomial is computed using the formula from
        Cordova-Gaiotto-Shao (Eq. C.28):

            Q̃_{[x],[y]} = Σ_{z ∈ [y]} Q_{x̄, z} · (-1)^{ℓ(x̄)} · (-1)^{ℓ(z)}

        where x̄ is the maximal representative of [x].

        Parameters
        ----------
        coset_x, coset_y : CosetRepresentative
            Coset representatives
        at_one : bool
            If True (default), evaluate at q=1

        Returns
        -------
        int or polynomial
            Q̃_{[x],[y]}(1) or Q̃_{[x],[y]}(q)

        Notes
        -----
        This formula is from Cordova, Gaiotto, Shao "Infrared Computations
        of Defect Schur Indices", Eq. (C.28).
        """
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self._W)
        parabolic = coset_x._parabolic

        # Get minimal representatives
        x_min = coset_x.representative
        y_min = coset_y.representative

        # Check Bruhat order on cosets
        if not bruhat.le(x_min, y_min):
            return 0

        # Get maximal representative of [x]
        x_max = self._maximal_coset_representative(x_min, parabolic)

        # Sum over elements in coset [y]
        result = 0
        for z in self._coset_elements(y_min, parabolic):
            if bruhat.le(x_max, z):
                q_val = self.Q_tilde(x_max, z, at_one=at_one)
                sign = (-1) ** (bruhat.length(x_max) + bruhat.length(z))
                result += sign * q_val

        return result

    def _maximal_coset_representative(self, w_min: Any, parabolic: "ParabolicSubgroup") -> Any:
        """Find the maximal length representative of a coset."""
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self._W)

        # Start from minimal rep and multiply by longest element of W_I
        # The maximal rep is w_min * w_I^0 where w_I^0 is longest in W_I
        current = w_min
        changed = True
        while changed:
            changed = False
            for i in parabolic.generators:
                s_i = self._W.simple_reflection(i)
                product = current * s_i
                if bruhat.length(product) > bruhat.length(current):
                    current = product
                    changed = True
                    break
        return current

    def _coset_elements(self, w_min: Any, parabolic: "ParabolicSubgroup") -> List[Any]:
        """Get all elements in the coset of w_min."""
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self._W)

        # Generate coset by multiplying w_min by all elements of W_I
        result = [w_min]
        queue = [w_min]
        visited = {self._element_key(w_min)}

        while queue:
            current = queue.pop(0)
            for i in parabolic.generators:
                s_i = self._W.simple_reflection(i)
                # Right multiplication for right cosets
                new_elem = current * s_i
                key = self._element_key(new_elem)
                if key not in visited:
                    visited.add(key)
                    result.append(new_elem)
                    queue.append(new_elem)

        return result

    # =========================================================================
    # Caching and Persistence
    # =========================================================================

    def save_cache(self, filename: Optional[str] = None) -> Path:
        """
        Save computed polynomials to disk.

        Parameters
        ----------
        filename : str, optional
            Filename (default: based on Cartan type)

        Returns
        -------
        Path
            Path to saved cache file
        """
        if filename is None:
            ct_str = str(self._cartan_type).replace(" ", "_")
            filename = f"kl_cache_{ct_str}.json"

        filepath = self._cache_dir / filename

        # Convert cache to serializable format
        cache_data = {
            "cartan_type": str(self._cartan_type),
            "Q_cache": {f"{k[0]}|{k[1]}": str(v) for k, v in self._Q_cache.items()},
            "P_cache": {f"{k[0]}|{k[1]}": str(v) for k, v in self._P_cache.items()},
        }

        with open(filepath, "w") as f:
            json.dump(cache_data, f, indent=2)

        return filepath

    def load_cache(self, filename: Optional[str] = None) -> bool:
        """
        Load cached polynomials from disk.

        Parameters
        ----------
        filename : str, optional
            Filename (default: based on Cartan type)

        Returns
        -------
        bool
            True if cache was loaded successfully
        """
        if filename is None:
            ct_str = str(self._cartan_type).replace(" ", "_")
            filename = f"kl_cache_{ct_str}.json"

        filepath = self._cache_dir / filename

        if not filepath.exists():
            return False

        try:
            with open(filepath) as f:
                cache_data = json.load(f)

            # Verify Cartan type matches
            if cache_data.get("cartan_type") != str(self._cartan_type):
                return False

            # Load Q cache
            for key_str, val_str in cache_data.get("Q_cache", {}).items():
                parts = key_str.split("|")
                if len(parts) == 2:
                    key = (
                        tuple(map(int, parts[0].split(",") if parts[0] else [])),
                        tuple(map(int, parts[1].split(",") if parts[1] else [])),
                    )
                    self._Q_cache[key] = SR(val_str)

            return True
        except Exception:
            return False

    # =========================================================================
    # Internal Methods
    # =========================================================================

    def _ensure_element(self, w: Any) -> Any:
        """Ensure w is an element of self._W."""
        if hasattr(w, "parent") and w.parent() == self._W:
            return w
        if isinstance(w, (list, tuple)):
            return self._W.from_reduced_word(w)
        return w

    def _element_key(self, w: Any) -> Tuple[int, ...]:
        """Get a hashable key for a Weyl group element."""
        if hasattr(w, "reduced_word"):
            return tuple(w.reduced_word())
        return (0,)  # Identity


@dataclass
class KLPolynomialCache:
    """
    Persistent cache for KL polynomial computations.

    This class manages disk-based caching of expensive KL polynomial
    computations, enabling reuse across sessions.

    Parameters
    ----------
    cache_dir : Path
        Directory for cache files
    cartan_type : str
        String representation of Cartan type
    """

    cache_dir: Path
    cartan_type: str
    _data: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Load existing cache if available."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._load()

    @property
    def filepath(self) -> Path:
        """Path to cache file."""
        safe_name = self.cartan_type.replace(" ", "_").replace("[", "").replace("]", "")
        return self.cache_dir / f"kl_{safe_name}.dat"

    def get(self, x_key: Tuple, y_key: Tuple) -> Optional[Any]:
        """Get cached value."""
        key = f"{x_key}|{y_key}"
        return self._data.get(key)

    def set(self, x_key: Tuple, y_key: Tuple, value: Any) -> None:
        """Set cached value."""
        key = f"{x_key}|{y_key}"
        self._data[key] = value

    def _load(self) -> None:
        """Load cache from disk."""
        if self.filepath.exists():
            try:
                with open(self.filepath) as f:
                    self._data = json.load(f)
            except Exception:
                self._data = {}

    def save(self) -> None:
        """Save cache to disk."""
        with open(self.filepath, "w") as f:
            json.dump(self._data, f)
