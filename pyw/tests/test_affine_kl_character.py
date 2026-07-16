import pytest


@pytest.mark.sage
def test_character_numerator_a21():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    legacy = kl_char.character_weight_list(lam, order=1)

    assert legacy
    assert all(len(entry) == 1 for entry in legacy)
    assert all(list(entry.values())[0] != 0 for entry in legacy)


@pytest.mark.sage
def test_character_numerator_accepts_explicit_translations():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)
    beta = ala.affine_weyl_group()._finite_coroot_space.simple_roots()[1]

    legacy = kl_char.character_weight_list(lam, order=1, translations=[0, beta])

    assert legacy
    assert all(list(entry.values())[0] != 0 for entry in legacy)


@pytest.mark.sage
def test_prepare_bounded_kl_orbit_does_not_evaluate_q_tilde(monkeypatch):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)
    character = KazhdanLusztigCharacter(algebra)

    def fail_q_tilde(*args, **kwargs):
        raise AssertionError("orbit preparation must not evaluate Q_tilde")

    monkeypatch.setattr(character.kl_polynomial, "Q_tilde", fail_q_tilde)
    orbit = character.prepare_bounded_kl_orbit(lambda_hat, order=1)

    assert orbit.representatives
    assert len(orbit.representatives) == len(orbit.weyl_to_weight)
    assert {
        tuple(int(index) for index in representative.reduced_word())
        for representative in orbit.representatives
    } == set(orbit.weyl_to_weight)


@pytest.mark.sage
def test_character_weight_list_consumes_prepared_orbit(monkeypatch):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)
    character = KazhdanLusztigCharacter(algebra)
    orbit = character.prepare_bounded_kl_orbit(lambda_hat, order=1)
    calls = []
    prepare_calls = []

    def prepared_orbit(*args, **kwargs):
        prepare_calls.append((args, kwargs))
        return orbit

    def constant_q_tilde(first, second, *, stabilizer_candidates):
        calls.append((first, second, tuple(stabilizer_candidates)))
        return 0 if len(calls) == 1 else 1

    monkeypatch.setattr(character, "prepare_bounded_kl_orbit", prepared_orbit)
    monkeypatch.setattr(character.kl_polynomial, "Q_tilde", constant_q_tilde)
    result = character.character_weight_list(
        lambda_hat,
        order=1,
        Lambda_hat=orbit.Lambda_hat,
        w_to_lambda=orbit.w_to_lambda,
        translations=("translation-sentinel",),
        show_progress=True,
        debug=True,
    )

    assert prepare_calls == [
        (
            (lambda_hat,),
            {
                "order": 1,
                "Lambda_hat": orbit.Lambda_hat,
                "w_to_lambda": orbit.w_to_lambda,
                "translations": ("translation-sentinel",),
                "show_progress": True,
                "debug": True,
            },
        )
    ]
    assert [call[1] for call in calls] == list(orbit.representatives)
    assert all(call[0] == orbit.w_to_lambda for call in calls)
    assert all(call[2] == orbit.stabilizer for call in calls)
    assert [next(iter(term)) for term in result] == [
        orbit.weyl_to_weight[tuple(int(index) for index in representative.reduced_word())]
        for representative in orbit.representatives[1:]
    ]


@pytest.mark.sage
def test_prepare_bounded_kl_orbit_requires_explicit_pair():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    algebra = AffineLieAlgebra(["A", 2, 1])
    lambda_hat = AffineWeight.affine_fundamental_weight(algebra, 1)
    character = KazhdanLusztigCharacter(algebra)

    with pytest.raises(ValueError, match="must be provided together"):
        character.prepare_bounded_kl_orbit(
            lambda_hat,
            order=0,
            Lambda_hat=lambda_hat,
        )


@pytest.mark.sage
def test_character_a21():
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.affine_weight import AffineWeight
    from pyw.core.character import KazhdanLusztigCharacter

    ala = AffineLieAlgebra(["A", 2, 1])
    lam = AffineWeight.affine_fundamental_weight(ala, 1)
    kl_char = KazhdanLusztigCharacter(ala)

    result = kl_char.character(lam, order=1)

    assert result != 0
