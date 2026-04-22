import contextlib
import importlib.util
import io
from pathlib import Path

import pytest


def _word_tuple(element):
    if hasattr(element, "affine_word"):
        return tuple(int(i) for i in element.affine_word())
    return tuple(int(i) for i in element.reduced_word())


def _load_legacy_algebra_module(repo_root: Path):
    legacy_path = repo_root / "demos" / "Algebra.py"
    spec = importlib.util.spec_from_file_location("legacy_demo_algebra", legacy_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.sage
@pytest.mark.parametrize(
    "cartan, coeff, index",
    [
        (["A", 2, 1], -2, 0),
        (["A", 2, 1], -2, 1),
        (["A", 2, 1], -1, 2),
        (["A", 2, 1], 2, 1),
        (["D", 4, 1], -2, 0),
        (["D", 4, 1], -2, 1),
        (["D", 4, 1], -1, 2),
        (["D", 4, 1], 2, 1),
    ],
)
def test_affine_kl_context_matches_legacy_translation_and_candidate_counts(cartan, coeff, index):
    from pyw.core.affine_lie_algebra import AffineLieAlgebra
    from pyw.core.character import KazhdanLusztigCharacter

    repo_root = Path(__file__).resolve().parents[2]
    legacy = _load_legacy_algebra_module(repo_root)

    ala = AffineLieAlgebra(cartan)
    kl_char = KazhdanLusztigCharacter(ala)
    fw = ala.fundamental_weights()
    lam = coeff * fw[index]

    with contextlib.redirect_stdout(io.StringIO()):
        legacy_alg = legacy.Alg(cartan, QLoad=False)
        legacy_lambda = coeff * legacy_alg.omega[index]
        legacy_Lambda = legacy_alg.GetLambda(legacy_lambda)
        legacy_selected_weight = legacy_Lambda + legacy_alg.rho
        order_min = int(
            legacy_alg.Tolambdakn(legacy_selected_weight)[-1]
            - legacy_alg.Tolambdakn(legacy_lambda)[-1]
        )

    for order in (0, 1, 2):
        with contextlib.redirect_stdout(io.StringIO()):
            legacy_translations = legacy_alg.GetTranslationsBynShift(
                legacy_selected_weight,
                order=order_min + order,
                order_min=None,
            )
            legacy_weyl = legacy_alg.GetWeylGroupForqSeries(order=order, T=legacy_translations)

        legacy_counts = (
            len({_word_tuple(t) for t in legacy_translations}),
            len({_word_tuple(w) for w in legacy_weyl}),
        )

        context = kl_char.build_context(lam, order=order)
        current_counts = (len(context.translations), len(context.weyl_lambda_candidates))

        assert current_counts == legacy_counts
