import pytest


@pytest.mark.sage
def test_bounded_orbit_cache_round_trip_reconstructs_exact_data(tmp_path, monkeypatch):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight, affine_weight_key
    from pyw.core.character import KazhdanLusztigCharacter
    from pyw.core.weyl_group import element_word_list

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)
    first_character = KazhdanLusztigCharacter(algebra)
    first_orbit = first_character.prepare_bounded_kl_orbit(
        lambda_hat, order=1, orbit_cache_dir=tmp_path
    )

    second_character = KazhdanLusztigCharacter(algebra)
    monkeypatch.setattr(
        second_character,
        "_translations_by_n_shift",
        lambda *_args, **_kwargs: pytest.fail("cache hit enumerated translations"),
    )
    second_orbit = second_character.prepare_bounded_kl_orbit(
        lambda_hat, order=1, orbit_cache_dir=tmp_path
    )

    assert [element_word_list(w) for w in second_orbit.representatives] == [
        element_word_list(w) for w in first_orbit.representatives
    ]
    assert [element_word_list(w) for w in second_orbit.stabilizer] == [
        element_word_list(w) for w in first_orbit.stabilizer
    ]
    assert {
        word: affine_weight_key(weight) for word, weight in second_orbit.weyl_to_weight.items()
    } == {
        word: affine_weight_key(weight) for word, weight in first_orbit.weyl_to_weight.items()
    }


@pytest.mark.sage
def test_bounded_orbit_cache_rejects_caller_provided_translations(tmp_path):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)

    with pytest.raises(ValueError, match="caller-provided translations"):
        KazhdanLusztigCharacter(algebra).prepare_bounded_kl_orbit(
            lambda_hat,
            order=1,
            translations=[],
            orbit_cache_dir=tmp_path,
        )
