from types import MappingProxyType

import pytest


def _small_orbit_and_partition():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight, affine_weight_key
    from pyw.core.character import BoundedKLOrbit
    from pyw.core.hybrid_affine_character import BoundedAffineVermaPartition
    from pyw.core.weyl_group import element_word_list

    algebra = AffineLieAlgebra(["A", 1, 1])
    highest_weight = AffineWeight.affine_fundamental_weight(algebra, 1)
    affine_weyl_group = algebra.affine_weyl_group_sage()
    representatives = (
        affine_weyl_group.one(),
        affine_weyl_group.from_reduced_word([0]),
    )
    simple_roots = algebra.affine_simple_roots()
    orbit_weights = (highest_weight, highest_weight + simple_roots[0])
    weyl_to_weight = {
        tuple(element_word_list(representative)): weight
        for representative, weight in zip(representatives, orbit_weights, strict=True)
    }
    orbit = BoundedKLOrbit(
        lambda_hat=highest_weight,
        Lambda_hat=highest_weight,
        w_to_lambda=affine_weyl_group.one(),
        stabilizer=(affine_weyl_group.one(),),
        representatives=representatives,
        weyl_to_weight=MappingProxyType(weyl_to_weight),
        weight_to_weyl=MappingProxyType(
            {
                affine_weight_key(weight): representative
                for representative, weight in zip(
                    representatives, orbit_weights, strict=True
                )
            }
        ),
        by_grade=MappingProxyType({}),
    )
    partition = BoundedAffineVermaPartition(algebra, (1, 1))
    return highest_weight, orbit, partition


@pytest.mark.sage
def test_kl_coefficient_filters_partition_and_reuses_q_tilde_cache():
    from pyw.core.hybrid_affine_character import (
        HybridCharacterStats,
        kl_coefficient_at_weight,
    )

    target, orbit, partition = _small_orbit_and_partition()

    class FakeKLPolynomial:
        def __init__(self):
            self.calls = []

        def Q_tilde(self, _x, y, *, stabilizer_candidates):  # noqa: N802
            self.calls.append((tuple(y.reduced_word()), tuple(stabilizer_candidates)))
            return 3

    kl_polynomial = FakeKLPolynomial()
    cache = {}
    stats = HybridCharacterStats()

    first = kl_coefficient_at_weight(
        target,
        orbit,
        partition,
        kl_polynomial,
        q_tilde_cache=cache,
        stats=stats,
    )
    second = kl_coefficient_at_weight(
        target,
        orbit,
        partition,
        kl_polynomial,
        q_tilde_cache=cache,
        stats=stats,
    )

    assert first == second == 6
    assert len(kl_polynomial.calls) == 2
    assert stats.q_tilde_requests == 4
    assert stats.q_tilde_unique_evaluations == 2
    assert stats.q_tilde_cache_hits == 2


@pytest.mark.sage
def test_hybrid_character_assembles_explicit_weight_coefficients(monkeypatch):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight, affine_weight_key
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter

    algebra = AffineLieAlgebra(["A", 1, 1])
    highest_weight = AffineWeight.affine_fundamental_weight(algebra, 1)

    class FakeKLCharacter:
        def __init__(self):
            self.algebra = algebra
            self.kl_polynomial = object()

        def prepare_bounded_kl_orbit(self, _weight, *, order):
            assert order == 1
            return object()

    engine = KazhdanLusztigFreudenthalCharacter(
        algebra, kl_character=FakeKLCharacter()
    )
    tree = __import__(
        "pyw.core.hybrid_affine_character", fromlist=["build_bounded_weight_tree"]
    ).build_bounded_weight_tree(highest_weight, order=1)
    dominant_weights = {
        affine_weight_key(weight): weight
        for weight in tree.weight_by_key.values()
        if weight.finite_part.is_dominant()
    }
    fake_multiplicities = {key: index + 1 for index, key in enumerate(dominant_weights)}

    monkeypatch.setattr(engine, "_degenerate_targets", lambda _tree: ())
    monkeypatch.setattr(
        "pyw.core.hybrid_affine_character.compute_bounded_freudenthal_multiplicities",
        lambda _tree, *, degenerate_multiplicity: fake_multiplicities,
    )

    character = engine.character(highest_weight, order=1)

    assert set(character) == {0, 1}
    assert sum(len(weights) for weights in character.values()) == sum(
        len(tree.weights_at_depth(depth)) for depth in range(2)
    )
    assert all(multiplicity > 0 for weights in character.values() for multiplicity in weights.values())


@pytest.mark.sage
@pytest.mark.slow
def test_d4_hybrid_character_uses_lazy_kl_and_matches_order_four_dimensions(tmp_path):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter
    from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

    algebra = AffineLieAlgebra(["D", 4, 1])
    highest_weight = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    kl_character = KazhdanLusztigCharacter(algebra)
    kl_character.kl_polynomial = KazhdanLusztigPolynomials(
        algebra.affine_weyl_group_sage(),
        cache_dir=tmp_path,
        persistent_cache=False,
    )
    engine = KazhdanLusztigFreudenthalCharacter(
        algebra, kl_character=kl_character
    )

    character = engine.character(highest_weight, order=4)
    dimensions = [sum(character[depth].values()) for depth in range(5)]
    stats = engine.profile_stats()

    assert dimensions == [1, 28, 329, 2632, 16380]
    assert stats["degenerate_targets"] == 7
    assert stats["candidate_query_hits"] == 238
    assert stats["q_tilde_requests"] == 238
    assert stats["q_tilde_unique_evaluations"] == 81
    assert stats["q_tilde_cache_hits"] == 157


@pytest.mark.sage
@pytest.mark.slow
def test_d4_hybrid_character_reuses_underlying_q_cache_in_process(tmp_path):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter
    from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

    algebra = AffineLieAlgebra(["D", 4, 1])
    highest_weight = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    kl_character = KazhdanLusztigCharacter(algebra)
    kl_character.kl_polynomial = KazhdanLusztigPolynomials(
        algebra.affine_weyl_group_sage(),
        cache_dir=tmp_path,
        persistent_cache=False,
    )
    engine = KazhdanLusztigFreudenthalCharacter(
        algebra, kl_character=kl_character
    )

    engine.character(highest_weight, order=4)
    engine.character(highest_weight, order=4)
    stats = engine.profile_stats()

    assert stats["kl_q_calls"] > 0
    assert stats["kl_q_cache_hits"] == stats["kl_q_calls"]
    assert stats["kl_q_tilde_calls"] == stats["q_tilde_unique_evaluations"] == 81
