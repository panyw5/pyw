"""Measure isolated cold and warm D4 hybrid-character calculations."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight
from pyw.core.character import KazhdanLusztigCharacter
from pyw.core.hybrid_affine_character import KazhdanLusztigFreudenthalCharacter
from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials


def run(cache_dir: Path, *, persistent_cache: bool, runs: int) -> list[dict[str, object]]:
    algebra = AffineLieAlgebra(["D", 4, 1])
    lambda_hat = -2 * AffineWeight.affine_fundamental_weight(algebra, 0)
    kl_character = KazhdanLusztigCharacter(algebra)
    kl_character.kl_polynomial = KazhdanLusztigPolynomials(
        algebra.affine_weyl_group_sage(),
        cache_dir=cache_dir,
        persistent_cache=persistent_cache,
    )
    engine = KazhdanLusztigFreudenthalCharacter(
        algebra,
        kl_character=kl_character,
        orbit_cache_dir=cache_dir / "orbit",
    )
    measurements = []
    for _ in range(runs):
        character = engine.character(lambda_hat, order=4)
        measurements.append(
            {
                "dimensions": [sum(character[depth].values()) for depth in range(5)],
                "stats": engine.profile_stats(),
            }
        )
    return measurements


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--persistent-cache", action="store_true")
    parser.add_argument("--runs", type=int, default=1)
    arguments = parser.parse_args()
    if arguments.cache_dir is None:
        with tempfile.TemporaryDirectory() as directory:
            print(
                json.dumps(
                    run(
                        Path(directory),
                        persistent_cache=arguments.persistent_cache,
                        runs=arguments.runs,
                    )
                )
            )
        return
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    print(
        json.dumps(
            run(
                arguments.cache_dir,
                persistent_cache=arguments.persistent_cache,
                runs=arguments.runs,
            )
        )
    )


if __name__ == "__main__":
    main()
