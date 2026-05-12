"""
Kazhdan-Lusztig Polynomials for Coxeter Groups

This module provides computation of Kazhdan-Lusztig polynomials and their
inverse/parabolic variants, essential for character formulas of admissible
representations.

Key features:
    - Standard KL polynomials P_{x,y}(q)
    - Inverse KL polynomials Q_{x,y}(q) on Weyl-group elements
    - Quotient/parabolic inverse KL polynomials Q̃_{[x],[y]}(q) on cosets
    - Integration with SageMath and coxeter3
    - Caching for expensive computations

The quotient inverse KL polynomial Q̃_{[x],[y]}(1) appears in the Kazhdan-Lusztig
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
import os
import shutil
import time
import warnings
from glob import glob
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Tuple, Union

from sage.all import QQ, RDF, SR, WeylGroup, matrix, var

if TYPE_CHECKING:
    from .bruhat import BruhatOrder, CosetRepresentative, ParabolicSubgroup
    from .affine_lie_algebra import AffineLieAlgebra


class KazhdanLusztigPolynomials:
    """
    Compute Kazhdan-Lusztig polynomials for Coxeter/Weyl groups.

    This class provides access to standard KL polynomials P_{x,y}(q),
    ordinary inverse KL polynomials Q_{x,y}(q) on Weyl-group elements,
    and quotient/parabolic inverse KL polynomials Q̃_{[x],[y]}(q) on cosets.

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
    The ordinary inverse KL polynomial Q_{x,y}(q) is related to P_{x,y}(q) by:
        Σ_z P_{x,z}(q) Q̃_{z,y}(q) = δ_{x,y}

    For the character formula on quotient data, we need Q̃_{[x],[y]}(1).
    """

    def __init__(
        self,
        coxeter_group: Any,
        cache_dir: Optional[Path] = None,
        persistent_cache: bool = True,
        auto_load_cache: bool = True,
        auto_save_cache: bool = True,
        auto_save_min_new_entries: int = 50,
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
        self.cartan_type = coxeter_group.cartan_type()

        try:
            from sage.all import CoxeterGroup as SageCoxeterGroup

            self.weyl_group = SageCoxeterGroup(self.cartan_type, implementation="coxeter3")
        except Exception:
            self.weyl_group = coxeter_group

        # Setup caching
        if cache_dir is None:
            cache_dir = Path.home() / ".pyw" / "kl_cache"
        self._cache_dir = cache_dir
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._persistent_cache_enabled = persistent_cache
        self._persistent_cache_auto_load = auto_load_cache
        self._persistent_cache_auto_save = auto_save_cache
        self._persistent_cache_auto_save_min_new_entries = max(1, int(auto_save_min_new_entries))
        self._persistent_cache_loaded = False
        self._persistent_cache_dirty = False
        self._persistent_cache_new_entries = 0

        # In-memory cache
        self._P_cache: Dict[Tuple[Any, Any], Any] = {}
        self._Q_cache: Dict[Tuple[Any, Any], Any] = {}
        self._invpol_cache: Dict[Tuple[Any, Any], Any] = {}
        self._Q_at_one_cache: Dict[Tuple[Any, Any], Any] = {}
        self._legacy_Q_cache: Dict[Tuple[Tuple[int, ...], Tuple[int, ...]], Any] = {}
        self._profiling_enabled = False
        self._profile_stats: Dict[str, Any] = {}

        # Initialize SageMath KL calculator
        self._sage_kl = None
        self._coxeter3 = None

        # Cython coxeter3 backend (direct C calls, no pexpect I/O)
        self._cython_group: Any = None
        self._cython_element_cache: Dict[tuple[int, ...], Any] = {}

        self._setup_backends()
        self._init_cython_backend()
        self.reset_profile_stats()
        if self._persistent_cache_enabled and self._persistent_cache_auto_load:
            self.load_cache()

    def reset_profile_stats(self) -> None:
        self._profile_stats = {
            "Q_calls": 0,
            "Q_total_seconds": 0.0,
            "Q_invpol_seconds": 0.0,
            "Q_invpol_calls": 0,
            "Q_cache_hits_at_one": 0,
            "Q_cache_hits_poly": 0,
            "Q_tilde_calls": 0,
            "Q_tilde_total_seconds": 0.0,
            "Q_tilde_stabilizer_terms": 0,
        }

    def set_profiling(self, enabled: bool) -> None:
        self._profiling_enabled = enabled

    def profile_stats(self) -> Dict[str, Any]:
        return dict(self._profile_stats)

    def _default_cache_filename(self) -> str:
        ct_str = str(self.cartan_type).replace(" ", "_")
        return f"kl_cache_{ct_str}.json"

    def _mark_persistent_cache_dirty(self, new_entries: int = 1) -> None:
        if not self._persistent_cache_enabled:
            return
        self._persistent_cache_dirty = True
        self._persistent_cache_new_entries += max(0, int(new_entries))
        self._maybe_autosave_cache()

    def _maybe_autosave_cache(self) -> None:
        if not self._persistent_cache_enabled:
            return
        if not self._persistent_cache_auto_save:
            return
        if self._persistent_cache_new_entries < self._persistent_cache_auto_save_min_new_entries:
            return
        self.save_cache()

    def _setup_backends(self) -> None:
        """Initialize computation backends."""
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        R = PolynomialRing(QQ, "q")
        self._q = R.gen()

        try:
            from sage.all import CoxeterGroup as SageCoxeterGroup

            ct = self.weyl_group.cartan_type()
            self._coxeter3_group = SageCoxeterGroup(ct, implementation="coxeter3")
        except Exception:
            self._coxeter3_group = None

        try:
            from sage.combinat.kazhdan_lusztig import KazhdanLusztigPolynomial

            self._sage_kl = KazhdanLusztigPolynomial(self.weyl_group, self._q)
        except ImportError:
            pass

        try:
            from coxeter3_sage import Coxeter3

            command = self._discover_coxeter_command()
            self._coxeter3 = Coxeter3(self.weyl_group, self._q, command=command)
        except (ImportError, Exception):
            pass

    def _init_cython_backend(self) -> None:
        """Initialize Cython coxeter3 backend for direct C calls (no pexpect I/O)."""
        try:
            from sage.libs.coxeter3.coxeter import get_CoxGroup as _CoxGroup, CoxGroupElement

            ct = self.cartan_type
            self._cython_group = _CoxGroup(ct)
            self._CoxGroupElement = CoxGroupElement
        except Exception:
            self._cython_group = None

    def _to_cython_element(self, w: Any) -> Any:
        """Convert a Sage Weyl element to a Cython CoxGroupElement (cached)."""
        if self._cython_group is None:
            return None
        key = self._word_tuple(w)
        try:
            return self._cython_element_cache[key]
        except KeyError:
            el = self._CoxGroupElement(self._cython_group, list(key))
            self._cython_element_cache[key] = el
            return el

    def _cython_Q_at_one(self, x: Any, y: Any) -> Any:
        """Compute Q(x,y) at q=1 via Cython coxeter3.

        Constructs the P-matrix evaluated at q=1 over the Bruhat interval [x,y]
        and inverts it to obtain Q(x,y,1).
        """
        cx = self._to_cython_element(x)
        cy = self._to_cython_element(y)
        if cx is None or cy is None:
            return None
        if not cx.bruhat_le(cy):
            return 0
        interval = self._cython_group.bruhat_interval(cx, cy)
        n = len(interval)
        if n == 1:
            return 1
        # P-matrix at q=1: P_mat[i,j] = (-1)^(l(j)-l(i)) * P_{i,j}(1) for i <= j
        P_mat = matrix(QQ, n, n)
        for i in range(n):
            li = len(interval[i])
            for j in range(i, n):
                if interval[i].bruhat_le(interval[j]):
                    lj = len(interval[j])
                    sign = 1 if (lj - li) % 2 == 0 else -1
                    p_val = self.P(
                        self.weyl_group.from_reduced_word(list(interval[i])),
                        self.weyl_group.from_reduced_word(list(interval[j])),
                        at_one=True,
                    )
                    P_mat[i, j] = sign * p_val
        Q_mat = P_mat.inverse()
        xw = self._word_tuple(x)
        yw = self._word_tuple(y)
        ix = next(i for i, w in enumerate(interval) if self._word_tuple(w) == xw)
        iy = next(i for i, w in enumerate(interval) if self._word_tuple(w) == yw)
        return Q_mat[ix, iy]

    def _discover_coxeter_command(self) -> str:
        env_command = os.environ.get("PYW_COXETER_COMMAND")
        if env_command:
            return env_command

        direct = shutil.which("coxeter")
        if direct:
            return direct

        candidates = sorted(glob("/private/var/tmp/sage-*/local/bin/coxeter"), reverse=True)
        for candidate in candidates:
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate

        return "coxeter"

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
        x = self._to_coxeter3(x)
        y = self._to_coxeter3(y)

        # Check cache
        cache_key = (self._element_key(x), self._element_key(y))
        if cache_key in self._P_cache:
            p = self._P_cache[cache_key]
            return p.subs({self._q: 1}) if at_one else p

        if self._coxeter3_group is not None:
            word = self._word_tuple(x)
            x3 = self._coxeter3_group.from_reduced_word(word)
            word2 = self._word_tuple(y)
            y3 = self._coxeter3_group.from_reduced_word(word2)
            p = self._coxeter3_group.kazhdan_lusztig_polynomial(x3, y3)
            self._P_cache[cache_key] = p
            return int(p.subs({self._q: 1})) if at_one else p

        if self._coxeter3 is not None:
            x_word = [i + 1 for i in self._word_tuple(x)]
            y_word = [i + 1 for i in self._word_tuple(y)]
            p = self._coxeter3.P(x_word, y_word)
            self._P_cache[cache_key] = p
            return int(p.subs({self._q: 1})) if at_one else p

        if self._sage_kl is not None:
            p = self._sage_kl.P(x, y)
            self._P_cache[cache_key] = p
            return p.subs({self._q: 1}) if at_one else p

        raise NotImplementedError(
            "KL polynomial computation requires coxeter3 or SageMath's KazhdanLusztigPolynomial"
        )

    def _P_by_words_experiment(self, x_word: tuple, y_word: tuple, at_one: bool = True) -> Any:
        cache_key = (x_word, y_word)
        if cache_key in self._P_cache:
            p = self._P_cache[cache_key]
            return int(p.subs({self._q: 1})) if at_one else p

        if self._coxeter3_group is not None:
            x3 = self._coxeter3_group.from_reduced_word(x_word)
            y3 = self._coxeter3_group.from_reduced_word(y_word)
            p = self._coxeter3_group.kazhdan_lusztig_polynomial(x3, y3)
            self._P_cache[cache_key] = p
            return int(p.subs({self._q: 1})) if at_one else p

        x_el = self.weyl_group.from_reduced_word(x_word)
        y_el = self.weyl_group.from_reduced_word(y_word)
        return self.P(x_el, y_el, at_one=at_one)

    # =========================================================================
    # Inverse KL Polynomials Q_{x,y}(q)
    # =========================================================================

    def Q(self, x: Any, y: Any, at_one: bool = False) -> Any:
        """
        Compute the ordinary inverse Kazhdan-Lusztig polynomial Q_{x,y}(q).

        Q_{x,y}(q) satisfies:
            Σ_z P_{x,z}(q) Q̃_{z,y}(q) = δ_{x,y}

        This is the Weyl-group-element-level inverse KL object. It is distinct
        from the quotient/parabolic quantity :meth:`Q_tilde`, whose arguments
        are coset representatives.

        Parameters
        ----------
        x, y : Weyl group elements
            Elements with x ≤ y in Bruhat order
        at_one : bool
            If True, return Q_{x,y}(1)

        Returns
        -------
        polynomial or int
            Q_{x,y}(q) or Q_{x,y}(1)

        Examples
        --------
        >>> W = WeylGroup(['A', 2])
        >>> kl = KazhdanLusztigPolynomials(W)
        >>> kl.Q(W.one(), W.long_element(), at_one=True)
        1
        """
        q_started = time.perf_counter() if self._profiling_enabled else None

        def _coerce_exact(value: Any) -> Any:
            if isinstance(value, float):
                return QQ(value)
            if hasattr(value, "is_numeric"):
                try:
                    if value.is_numeric():
                        return QQ(value)
                except Exception:
                    pass
            if hasattr(value, "is_integer"):
                try:
                    if value.is_integer():
                        return int(value)
                except Exception:
                    pass
            if not hasattr(value, "polynomial"):
                return value
            try:
                poly_rdf = value.polynomial(RDF)
                poly_qq = QQ[poly_rdf.parent().variable_name()](poly_rdf)
                return SR(poly_qq)
            except Exception:
                return value

        # Check in-memory cache first (works for both Cython and pexpect paths).
        cache_key = (self._element_key(x), self._element_key(y))
        if at_one and cache_key in self._Q_at_one_cache:
            cached_value = _coerce_exact(self._Q_at_one_cache[cache_key])
            self._Q_at_one_cache[cache_key] = cached_value
            if self._profiling_enabled:
                self._profile_stats["Q_calls"] += 1
                self._profile_stats["Q_cache_hits_at_one"] += 1
                if q_started is not None:
                    self._profile_stats["Q_total_seconds"] += time.perf_counter() - q_started
            return cached_value
        if not at_one and cache_key in self._Q_cache:
            cached_value = _coerce_exact(self._Q_cache[cache_key])
            self._Q_cache[cache_key] = cached_value
            if self._profiling_enabled:
                self._profile_stats["Q_calls"] += 1
                self._profile_stats["Q_cache_hits_poly"] += 1
                if q_started is not None:
                    self._profile_stats["Q_total_seconds"] += time.perf_counter() - q_started
            return cached_value

        # Fast path: Cython coxeter3 (direct C, no pexpect I/O).
        if at_one and self._cython_group is not None:
            try:
                result = self._cython_Q_at_one(x, y)
                if result is not None:
                    result = _coerce_exact(result)
                    cache_key = (self._element_key(x), self._element_key(y))
                    self._Q_at_one_cache[cache_key] = result
                    self._mark_persistent_cache_dirty()
                    if self._profiling_enabled:
                        self._profile_stats["Q_calls"] += 1
                        if q_started is not None:
                            self._profile_stats["Q_total_seconds"] += (
                                time.perf_counter() - q_started
                            )
                    return result
            except Exception:
                pass  # fall through to existing backends

        x = self._to_coxeter3(x)
        y = self._to_coxeter3(y)

        if self._coxeter3 is not None and hasattr(self._coxeter3, "invpol"):
            try:
                invpol_started = time.perf_counter() if self._profiling_enabled else None
                result = self._coxeter3.invpol(x, y)
                result = _coerce_exact(result)
                if self._profiling_enabled:
                    self._profile_stats["Q_invpol_calls"] += 1
                    if invpol_started is not None:
                        self._profile_stats["Q_invpol_seconds"] += (
                            time.perf_counter() - invpol_started
                        )
                self._Q_cache[cache_key] = result
                if at_one:
                    value_at_one = result.subs({self._q: 1}) if hasattr(result, "subs") else result
                    value_at_one = _coerce_exact(value_at_one)
                    if cache_key not in self._Q_at_one_cache:
                        self._Q_at_one_cache[cache_key] = value_at_one
                        self._mark_persistent_cache_dirty()
                    else:
                        self._Q_at_one_cache[cache_key] = value_at_one
                    if self._profiling_enabled:
                        self._profile_stats["Q_calls"] += 1
                        if q_started is not None:
                            self._profile_stats["Q_total_seconds"] += (
                                time.perf_counter() - q_started
                            )
                    return value_at_one
                if self._profiling_enabled:
                    self._profile_stats["Q_calls"] += 1
                    if q_started is not None:
                        self._profile_stats["Q_total_seconds"] += time.perf_counter() - q_started
                return result
            except Exception:
                pass

        if not at_one:
            result = self._compute_inverse_kl_by_matrix_inversion(x, y, at_one=False)
            result = _coerce_exact(result)
            self._Q_cache[cache_key] = result
            if self._profiling_enabled:
                self._profile_stats["Q_calls"] += 1
                if q_started is not None:
                    self._profile_stats["Q_total_seconds"] += time.perf_counter() - q_started
            return result

        result = self._compute_inverse_kl_by_matrix_inversion(x, y, at_one=True)
        result = _coerce_exact(result)
        if cache_key not in self._Q_at_one_cache:
            self._Q_at_one_cache[cache_key] = result
            self._mark_persistent_cache_dirty()
        else:
            self._Q_at_one_cache[cache_key] = result
        if self._profiling_enabled:
            self._profile_stats["Q_calls"] += 1
            if q_started is not None:
                self._profile_stats["Q_total_seconds"] += time.perf_counter() - q_started
        return result

    def Q_tilde_experiment(
        self,
        coset_x: "CosetRepresentative",
        coset_y: "CosetRepresentative",
        at_one: bool = True,
    ) -> Any:
        """
        Compute the quotient/parabolic inverse KL polynomial Q̃_{[x],[y]}(q).

        This is the coset-level object attached to the quotient W / W_I. It is
        implemented by :meth:`parabolic_Q_tilde_experiment`, which currently supports
        right cosets in finite Weyl groups.
        """
        return self.parabolic_Q_tilde_experiment(coset_x, coset_y, at_one=at_one)

    def invpol_experiment(self, x: Any, y: Any) -> Any:
        """
        Compute inverse KL polynomial using coxeter3's invpol.

        This is a direct interface to coxeter3's invpol command,
        which computes the ordinary inverse KL polynomial Q_{x,y}(q)
        efficiently.

        Parameters
        ----------
        x, y : Weyl group elements
            Elements with x ≤ y in Bruhat order

        Returns
        -------
        polynomial
            Q_{x,y}(q) as a symbolic expression

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

        x = self._to_coxeter3(x)
        y = self._to_coxeter3(y)

        # Check cache
        cache_key = (self._element_key(x), self._element_key(y))
        if cache_key in self._invpol_cache:
            return self._invpol_cache[cache_key]

        # Call coxeter3
        result = self._coxeter3.invpol(x, y)
        self._invpol_cache[cache_key] = result
        return result

    def _compute_inverse_kl_by_matrix_inversion(self, x: Any, y: Any, at_one: bool) -> Any:
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self.weyl_group)

        if not bruhat.le(x, y):
            return 0

        interval = bruhat.interval(x, y)
        n = len(interval)

        if n == 1:
            return 1

        ring = QQ if at_one else self._q.parent()
        P_mat = matrix(ring, n, n)
        for i, w_i in enumerate(interval):
            for j, w_j in enumerate(interval):
                if bruhat.le(w_i, w_j):
                    P_mat[i, j] = self.P(w_i, w_j, at_one=at_one)

        try:
            Q_mat = P_mat.inverse()
        except Exception as e:
            raise RuntimeError(f"Failed to invert P-matrix: {e}")

        key_x = self._element_key(x)
        key_y = self._element_key(y)
        x_idx = next(i for i, w in enumerate(interval) if self._element_key(w) == key_x)
        y_idx = next(i for i, w in enumerate(interval) if self._element_key(w) == key_y)

        return Q_mat[x_idx, y_idx]

    # =========================================================================
    # Parabolic (Coset) KL Polynomials
    # =========================================================================

    def parabolic_Q_tilde_experiment(
        self,
        coset_x: "CosetRepresentative",
        coset_y: "CosetRepresentative",
        at_one: bool = False,
    ) -> Any:
        """
        Compute parabolic inverse KL polynomial Q̃_{[x],[y]}(q).

        For right cosets [x], [y] in W / W_I, the parabolic
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

        bruhat = BruhatOrder(self.weyl_group)
        parabolic = coset_x._parabolic

        if parabolic != coset_y._parabolic:
            raise ValueError("Coset inputs must use the same parabolic subgroup")
        if coset_x._left != coset_y._left:
            raise ValueError("Coset inputs must use the same left/right convention")
        if coset_x._left:
            raise NotImplementedError(
                "parabolic_Q_tilde_experiment currently supports right cosets only"
            )
        if not self.weyl_group.is_finite():
            raise ValueError("parabolic_Q_tilde_experiment currently supports finite groups only")

        # Get minimal representatives
        x_min = coset_x.representative
        y_min = coset_y.representative

        # Check Bruhat order on cosets
        if not bruhat.le(x_min, y_min):
            return 0

        # Get maximal representative of [x]
        x_max = self._maximal_representative_in_coset_experiment(x_min, parabolic)

        # Sum over elements in coset [y]
        result = 0
        for z in self._enumerate_coset_elements_experiment(y_min, parabolic):
            if bruhat.le(x_max, z):
                q_val = self.Q(x_max, z, at_one=at_one)
                sign = (-1) ** (bruhat.length(x_max) + bruhat.length(z))
                result += sign * q_val

        return result

    @staticmethod
    def _is_subseq(a: tuple, b: tuple) -> bool:
        """Return True if tuple *a* is a subsequence of tuple *b* (subword criterion)."""
        it = iter(b)
        return all(x in it for x in a)

    @staticmethod
    def _word_tuple(w: Any) -> tuple[int, ...]:
        if hasattr(w, "reduced_word"):
            try:
                return tuple(int(i) for i in w.reduced_word())
            except Exception:
                pass
        if hasattr(w, "reduced_word_list"):
            try:
                return tuple(int(i) for i in w.reduced_word_list())
            except Exception:
                pass
        try:
            return tuple(int(i) for i in w)
        except Exception:
            return ()

    @classmethod
    def _word_length(cls, w: Any) -> int:
        if hasattr(w, "length"):
            try:
                return int(w.length())
            except Exception:
                pass
        return len(cls._word_tuple(w))

    @classmethod
    def _build_word_cache(cls, elements: Iterable[Any]) -> Dict[int, tuple]:
        """Pre-compute reduced words for all elements, keyed by id()."""
        cache: Dict[int, tuple] = {}
        for e in elements:
            eid = id(e)
            if eid not in cache:
                cache[eid] = cls._word_tuple(e)
        return cache

    def _bruhat_le_by_words(
        self,
        x: Any,
        y: Any,
        word_cache: Optional[Dict[int, tuple]] = None,
    ) -> bool:
        """Bruhat order x ≤ y via subword criterion on reduced words.

        Significantly faster than SageMath's matrix-based ``bruhat_le`` for
        affine Weyl group elements because it avoids GAP matrix hashing.
        Falls back to SageMath ``bruhat_le`` when *word_cache* is not supplied.
        """
        if word_cache is None:
            return bool(x.bruhat_le(y))
        wx_cached = word_cache.get(id(x))
        wy_cached = word_cache.get(id(y))
        wx = wx_cached if wx_cached is not None else self._word_tuple(x)
        wy = wy_cached if wy_cached is not None else self._word_tuple(y)
        return self._is_subseq(wx, wy)

    def affine_bounded_interval_experiment(
        self,
        x: Any,
        y: Any,
        *,
        candidates: Iterable[Any],
        word_cache: Optional[Dict[int, tuple]] = None,
    ) -> List[Any]:
        """Bruhat interval restricted to an explicit bounded affine candidate set.

        Parameters
        ----------
        word_cache:
            Optional pre-computed ``{id(elem): reduced_word_tuple}`` mapping.
            When supplied, Bruhat comparisons use the fast subword criterion
            instead of SageMath's matrix-based ``bruhat_le``.
        """
        x_internal = self._to_coxeter3(x)
        y_internal = self._to_coxeter3(y)

        candidate_pairs = [(w, self._to_coxeter3(w)) for w in candidates]
        candidates_list = [internal for _, internal in candidate_pairs]
        if word_cache is None:
            word_cache = self._build_word_cache([x_internal, y_internal] + candidates_list)

        wx_cached = word_cache.get(id(x_internal))
        wy_cached = word_cache.get(id(y_internal))
        lx = len(wx_cached) if wx_cached is not None else self._word_length(x_internal)
        ly = len(wy_cached) if wy_cached is not None else self._word_length(y_internal)

        filtered_pairs = [
            (original, internal)
            for original, internal in candidate_pairs
            if lx <= len(word_cache.get(id(internal), ())) <= ly
            and self._bruhat_le_by_words(x_internal, internal, word_cache)
            and self._bruhat_le_by_words(internal, y_internal, word_cache)
        ]

        key_x = self._element_key(x_internal)
        key_y = self._element_key(y_internal)
        filtered_keys = {self._element_key(w) for _, w in filtered_pairs}
        if key_x not in filtered_keys or key_y not in filtered_keys:
            raise ValueError("bounded candidate set must contain both interval endpoints")

        return [original for original, _ in filtered_pairs]

    def affine_bounded_Q_experiment(
        self,
        x: Any,
        y: Any,
        *,
        candidates: Iterable[Any],
        at_one: bool = True,
        word_cache: Optional[Dict[int, tuple]] = None,
    ) -> Any:
        """Compute ordinary inverse KL polynomial on a bounded affine interval.

        For affine groups, direct interval enumeration is infinite globally but
        locally finite. This helper lets callers provide an explicit bounded
        candidate set coming from translation bounds, then performs the same
        matrix-inversion construction used by the finite fallback.

        Parameters
        ----------
        word_cache:
            Optional pre-computed reduced-word cache (see
            :meth:`affine_bounded_interval_experiment`).  Pass the same cache across all
            calls within one ``numerator_terms`` computation to avoid
            recomputing reduced words for every interval.
        """
        x = self._to_coxeter3(x)
        y = self._to_coxeter3(y)

        if word_cache is None:
            candidates_list = [self._to_coxeter3(w) for w in candidates]
            word_cache = self._build_word_cache([x, y] + candidates_list)

        if not self._bruhat_le_by_words(x, y, word_cache):
            return 0

        if self._coxeter3 is not None and hasattr(self._coxeter3, "invpol"):
            try:
                result = self._coxeter3.invpol(x, y)
                if at_one:
                    return result.subs({self._q: 1}) if hasattr(result, "subs") else result
                return result
            except Exception:
                pass

        candidates_list = [self._to_coxeter3(w) for w in candidates]
        interval = self.affine_bounded_interval_experiment(
            x, y, candidates=candidates_list, word_cache=word_cache
        )

        n = len(interval)
        if n == 1:
            return 1

        ring = QQ if at_one else self._q.parent()
        p_matrix = matrix(ring, n, n)
        for i, w_i in enumerate(interval):
            for j, w_j in enumerate(interval):
                if self._bruhat_le_by_words(w_i, w_j, word_cache):
                    xi = word_cache.get(id(w_i))
                    xj = word_cache.get(id(w_j))
                    if xi is not None and xj is not None:
                        p_matrix[i, j] = self._P_by_words_experiment(xi, xj, at_one=at_one)
                    else:
                        p_matrix[i, j] = self.P(w_i, w_j, at_one=at_one)

        q_matrix = p_matrix.inverse()
        key_x = self._element_key(x)
        key_y = self._element_key(y)
        idx_x = next(k for k, w in enumerate(interval) if self._element_key(w) == key_x)
        idx_y = next(k for k, w in enumerate(interval) if self._element_key(w) == key_y)
        return q_matrix[idx_x, idx_y]

    def affine_stabilizer_experiment(
        self,
        Lambda: Any,
        *,
        rho_hat: Any,
        candidates: Iterable[Any],
        algebra: Optional["AffineLieAlgebra"] = None,
    ) -> List[Any]:
        """Return the bounded affine stabilizer of ``Lambda`` under the dot action.

        When ``algebra`` is provided, candidate elements are interpreted through
        the semidirect-product affine Weyl wrapper so their action is evaluated
        on affine weights rather than Sage's default root-domain action.
        """
        from .affine_weight import AffineWeight

        result: List[Any] = []
        candidates_list = list(candidates)
        semidirect = algebra.affine_weyl_group() if algebra is not None else None
        target = Lambda + rho_hat
        target_affine = None
        rho_hat_affine = None
        Lambda_affine = None

        if algebra is not None:
            target_affine = AffineWeight.from_sagemath(algebra, target)
            rho_hat_affine = AffineWeight.from_sagemath(algebra, rho_hat)
            Lambda_affine = AffineWeight.from_sagemath(algebra, Lambda)

        for w in candidates_list:
            if semidirect is not None:
                if hasattr(w, "reduced_word"):
                    affine_word = w.word_list() if hasattr(w, "word_list") else self._word_tuple(w)
                    element = semidirect.from_word(affine_word)
                    stabilized = element.action(target_affine)
                else:
                    element = w
                    stabilized = element.action(target_affine)

                if stabilized - rho_hat_affine == Lambda_affine:
                    result.append(w)
                continue

            element = self._to_coxeter3(w)
            if element.action(target) - rho_hat == Lambda:
                result.append(w)

        identity = self.weyl_group.one()
        if candidates_list and hasattr(candidates_list[0], "parent"):
            try:
                identity = candidates_list[0].parent().one()
            except Exception:
                pass
        if all(self._element_key(w) != self._element_key(identity) for w in result):
            result.append(identity)

        return sorted(result, key=lambda w: (self._word_length(w), self._word_tuple(w)))

    def affine_bounded_parabolic_Q_tilde_experiment(
        self,
        x_min: Any,
        y_min: Any,
        *,
        candidates: Iterable[Any],
        stabilizer_candidates: Iterable[Any],
        at_one: bool = True,
        word_cache: Optional[Dict[int, tuple]] = None,
    ) -> Any:
        """Compute affine coset-level Q̃ using the direct legacy summation path.

        Notes
        -----
        Deprecated in favor of :meth:`Q_tilde` for default affine
        KL character assembly. Keep this wrapper only for compatibility with
        older call sites that still pass ``candidates``.

        ``candidates`` is accepted for API compatibility with older bounded
        affine workflows. It is intentionally not used to filter the coset
        summation domain.
        During benchmarking on the current ``KazhdanLusztigCharacter`` path,
        the previous bounded-candidates implementation was about 50x slower
        than the direct legacy summation because it repeatedly rebuilt bounded
        intervals and inverted their P-matrices.

        More importantly, filtering ``y_min * stabilizer`` through a bounded
        candidate set changes the mathematical object being computed: it turns
        a full coset sum into a candidate-truncated sum.  The direct path below
        matches the original ``_legacy_qtilde_at_one_direct`` behavior used by
        character assembly and should therefore be treated as the standard
        affine implementation unless a caller explicitly wants a truncated
        auxiliary computation. The only remaining use of ``candidates`` here is
        as a fallback support set for ordinary ``Q`` computations when the
        direct inverse-KL backend is unavailable.
        """
        warnings.warn(
            "affine_bounded_parabolic_Q_tilde_experiment() is deprecated for default KL character assembly; "
            "use Q_tilde() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        x_min = self._to_coxeter3(x_min, word_cache=word_cache)
        y_min = self._to_coxeter3(y_min, word_cache=word_cache)
        bounded_candidates = [self._to_coxeter3(w, word_cache=word_cache) for w in candidates]
        stabilizer = [self._to_coxeter3(w, word_cache=word_cache) for w in stabilizer_candidates]

        if word_cache is None:
            word_cache = self._build_word_cache([x_min, y_min] + bounded_candidates + stabilizer)

        if not self._bruhat_le_by_words(x_min, y_min, word_cache):
            return 0

        x_max = self._bounded_maximal_representative(x_min, stabilizer=stabilizer)
        if id(x_max) not in word_cache:
            word_cache[id(x_max)] = self._word_tuple(x_max)

        result = 0
        for stabilizer_element in stabilizer:
            coset_element = y_min * stabilizer_element
            if id(coset_element) not in word_cache:
                word_cache[id(coset_element)] = self._word_tuple(coset_element)

            if self._coxeter3 is not None and hasattr(self._coxeter3, "invpol"):
                q_value = self.Q(x_max, coset_element, at_one=at_one)
            else:
                ordinary_q_candidates_by_key = {
                    self._element_key(candidate): candidate
                    for candidate in bounded_candidates + [x_max, coset_element]
                }
                q_value = self.affine_bounded_Q_experiment(
                    x_max,
                    coset_element,
                    candidates=ordinary_q_candidates_by_key.values(),
                    at_one=at_one,
                    word_cache=word_cache,
                )

            sign = (-1) ** (len(word_cache[id(x_max)]) - len(word_cache[id(coset_element)]))
            result += sign * q_value

        if hasattr(result, "full_simplify"):
            result = result.full_simplify()
        return result

    def Q_tilde(
        self,
        x_min: Any,
        y_min: Any,
        *,
        stabilizer_candidates: Iterable[Any],
        at_one: bool = True,
        word_cache: Optional[Dict[int, tuple]] = None,
    ) -> Any:
        q_tilde_started = time.perf_counter() if self._profiling_enabled else None
        x_min = self._to_coxeter3(x_min, word_cache=word_cache)
        y_min = self._to_coxeter3(y_min, word_cache=word_cache)
        stabilizer = [self._to_coxeter3(w, word_cache=word_cache) for w in stabilizer_candidates]

        if word_cache is None:
            word_cache = self._build_word_cache([x_min, y_min] + stabilizer)

        # NOTE: Legacy Q_tilde follows the MyAlgebra.py convention:
        # do NOT short-circuit when x_min not <= y_min in Bruhat.
        # The stabilizer sum over y_min * stabilizer can still yield
        # non-zero contributions (critical for matching legacy results).
        # Old code: Qtilde(x,y,subgroup) = sum_{s in subgroup}
        #   Q(MaxRep(x,subgroup), y*s) * (-1)^(l(xbar)-l(y*s))
        # with no pre-filter on x <= y.

        x_max = self._bounded_maximal_representative(x_min, stabilizer=stabilizer)
        if id(x_max) not in word_cache:
            word_cache[id(x_max)] = self._word_tuple(x_max)

        result = 0
        for stabilizer_element in stabilizer:
            coset_element = y_min * stabilizer_element
            if id(coset_element) not in word_cache:
                word_cache[id(coset_element)] = self._word_tuple(coset_element)
            q_value = self.Q(x_max, coset_element, at_one=at_one)
            sign = (-1) ** (len(word_cache[id(x_max)]) - len(word_cache[id(coset_element)]))
            result += sign * q_value

        if hasattr(result, "full_simplify"):
            result = result.full_simplify()
        if hasattr(result, "polynomial"):
            try:
                poly_rdf = result.polynomial(RDF)
                poly_qq = QQ[poly_rdf.parent().variable_name()](poly_rdf)
                result = SR(poly_qq)
            except Exception:
                pass
        if isinstance(result, float):
            result = QQ(result)
        elif hasattr(result, "is_numeric"):
            try:
                if result.is_numeric():
                    result = QQ(result)
            except Exception:
                pass
        if self._profiling_enabled:
            self._profile_stats["Q_tilde_calls"] += 1
            self._profile_stats["Q_tilde_stabilizer_terms"] += len(stabilizer)
            if q_tilde_started is not None:
                self._profile_stats["Q_tilde_total_seconds"] += (
                    time.perf_counter() - q_tilde_started
                )
        return result

    def _bounded_maximal_representative(self, w_min: Any, *, stabilizer: Iterable[Any]) -> Any:
        """Find the maximal representative inside a bounded right coset."""
        current = self._to_coxeter3(w_min)
        current_length = self._word_length(current)
        for s in stabilizer:
            candidate = current * self._to_coxeter3(s)
            candidate_length = self._word_length(candidate)
            if candidate_length > current_length:
                current = candidate
                current_length = candidate_length
        return current

    def _bounded_right_coset_elements_experiment(
        self,
        w_min: Any,
        *,
        stabilizer: Iterable[Any],
        candidates: Iterable[Any],
        candidate_set: Optional[Dict[Tuple[int, ...], Any]] = None,
        word_cache: Optional[Dict[int, tuple]] = None,
    ) -> List[Any]:
        """Enumerate right-coset elements present in a bounded candidate set."""
        if candidate_set is None:
            candidate_set = {
                self._element_key(self._to_coxeter3(w)): self._to_coxeter3(w)
                for w in candidates
            }
        result: Dict[Tuple[int, ...], Any] = {}
        base = self._to_coxeter3(w_min)
        for s in stabilizer:
            candidate = base * self._to_coxeter3(s)
            key = self._element_key(candidate)
            if key in candidate_set:
                result[key] = candidate_set[key]
        base_key = self._element_key(base)
        if base_key not in result and base_key in candidate_set:
            result[base_key] = candidate_set[base_key]
        if word_cache is not None:
            return sorted(result.values(), key=lambda w: len(word_cache.get(id(w), ())))
        return sorted(result.values(), key=lambda w: (self._word_length(w), self._word_tuple(w)))

    def _maximal_representative_in_coset_experiment(
        self, w_min: Any, parabolic: "ParabolicSubgroup"
    ) -> Any:
        """Find the maximal length representative of a coset."""
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self.weyl_group)

        # Start from minimal rep and multiply by longest element of W_I
        # The maximal rep is w_min * w_I^0 where w_I^0 is longest in W_I
        current = w_min
        changed = True
        while changed:
            changed = False
            for i in parabolic.generators:
                s_i = self.weyl_group.simple_reflection(i)
                product = current * s_i
                if bruhat.length(product) > bruhat.length(current):
                    current = product
                    changed = True
                    break
        return current

    def _enumerate_coset_elements_experiment(
        self, w_min: Any, parabolic: "ParabolicSubgroup"
    ) -> List[Any]:
        """Get all elements in the coset of w_min."""
        from .bruhat import BruhatOrder

        bruhat = BruhatOrder(self.weyl_group)

        # Generate coset by multiplying w_min by all elements of W_I
        result = [w_min]
        queue = [w_min]
        visited = {self._element_key(w_min)}

        while queue:
            current = queue.pop(0)
            for i in parabolic.generators:
                s_i = self.weyl_group.simple_reflection(i)
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

    def save_cache(self, filename: Optional[str] = None, force: bool = False) -> Path:
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
            filename = self._default_cache_filename()

        if not force and not self._persistent_cache_dirty and self._persistent_cache_loaded:
            return self._cache_dir / filename

        filepath = self._cache_dir / filename

        # Convert cache to serializable format
        cache_data = {
            "cache_version": self.CACHE_VERSION_experiment,
            "cartan_type": str(self.cartan_type),
            "value_kind": "Q_at_one",
            "Q_at_one_cache": {
                self._cache_key_to_string_experiment(k): self._json_scalar_experiment(v)
                for k, v in self._Q_at_one_cache.items()
            },
        }

        with open(filepath, "w") as f:
            json.dump(cache_data, f, indent=2)

        self._persistent_cache_dirty = False
        self._persistent_cache_new_entries = 0

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
            filename = self._default_cache_filename()

        filepath = self._cache_dir / filename

        if not filepath.exists():
            return False

        try:
            with open(filepath) as f:
                cache_data = json.load(f)

            if cache_data.get("cache_version") != self.CACHE_VERSION_experiment:
                return False

            if cache_data.get("cartan_type") != str(self.cartan_type):
                return False

            if cache_data.get("value_kind") not in {"Q_at_one", "Q_tilde_at_one"}:
                return False

            for key_str, value in cache_data.get("Q_at_one_cache", {}).items():
                self._Q_at_one_cache[self._parse_cache_key_experiment(key_str)] = (
                    self._parse_json_scalar_experiment(value)
                )

            self._persistent_cache_loaded = True
            self._persistent_cache_dirty = False
            self._persistent_cache_new_entries = 0
            return True
        except Exception:
            return False

    def save_cache_experiment(self, filename: Optional[str] = None) -> Path:
        return self.save_cache(filename=filename, force=True)

    def load_cache_experiment(self, filename: Optional[str] = None) -> bool:
        return self.load_cache(filename=filename)

    # =========================================================================
    # Internal Methods
    # =========================================================================

    def _to_coxeter3(self, w: Any, word_cache: Optional[Dict[int, tuple]] = None) -> Any:
        if hasattr(w, "parent") and w.parent() == self.weyl_group:
            return w
        if isinstance(w, (list, tuple)):
            return self.weyl_group.from_reduced_word(list(w))
        if word_cache is not None and id(w) in word_cache:
            return self.weyl_group.from_reduced_word(list(word_cache[id(w)]))
        word = self._word_tuple(w)
        if hasattr(w, "reduced_word") or hasattr(w, "reduced_word_list"):
            return self.weyl_group.from_reduced_word(list(word))
        if hasattr(w, "parent_group"):
            return self.weyl_group.from_reduced_word(list(word))
        return w

    def _element_key(self, w: Any) -> Tuple[int, ...]:
        """Get a hashable key for a Weyl group element."""
        return self._word_tuple(w)

    def _cache_key_to_string_experiment(
        self, cache_key: Tuple[Tuple[int, ...], Tuple[int, ...]]
    ) -> str:
        return json.dumps([list(cache_key[0]), list(cache_key[1])])

    def _parse_cache_key_experiment(self, key_str: str) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
        left, right = json.loads(key_str)
        return (tuple(left), tuple(right))

    def _json_scalar_experiment(self, value: Any) -> Any:
        if isinstance(value, bool) or value is None:
            return value
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            return int(value) if value.is_integer() else str(value)
        if isinstance(value, str):
            return value
        if hasattr(value, "numerator") and hasattr(value, "denominator"):
            try:
                numerator = (
                    int(value.numerator()) if callable(value.numerator) else int(value.numerator)
                )
                denominator = (
                    int(value.denominator())
                    if callable(value.denominator)
                    else int(value.denominator)
                )
                if denominator != 0:
                    return {"__qq__": [numerator, denominator]}
            except Exception:
                pass
        if hasattr(value, "is_integer") and value.is_integer():
            return int(value)
        try:
            if hasattr(value, "is_integer") and not value.is_integer():
                raise ValueError
            return int(value)
        except (TypeError, ValueError):
            return str(value)

    def _parse_json_scalar_experiment(self, value: Any) -> Any:
        if isinstance(value, dict) and "__qq__" in value:
            payload = value["__qq__"]
            if isinstance(payload, list) and len(payload) == 2:
                return QQ(payload[0]) / QQ(payload[1])
        if isinstance(value, float):
            return int(value) if value.is_integer() else QQ(value)
        return value

    CACHE_VERSION_experiment = 1
