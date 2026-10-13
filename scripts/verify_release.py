"""Verify distribution contents and an installed pyw in a SageMath environment.

Run with Python -I outside the source checkout after installing the built wheel.
This deliberately validates the minimum release scope, not the whole test suite.
"""

import argparse
import ast
import re
import shutil
import tarfile
import tempfile
import tomllib
from pathlib import Path
from zipfile import ZipFile


TEST_FILES = (
    "test_root_system.py",
    "test_coweight_conversion.py",
    "test_fundamental_coweights.py",
    "test_level.py",
    "test_admissible.py",
    "test_weight_space.py",
)
SUBPACKAGES = ("core", "fractional", "embedding", "algorithms", "utils")


def verify(project: Path, dist: Path, run_tests: bool) -> None:
    config = tomllib.loads((project / "pyproject.toml").read_text())
    version = config["project"]["version"]
    distribution_name = config["project"]["name"]
    init = ast.parse((project / "pyw" / "__init__.py").read_text())
    source_version = next(
        ast.literal_eval(node.value)
        for node in init.body
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "__version__" for t in node.targets)
    )
    assert version == source_version, "Source and metadata versions differ"
    wheels = list(dist.glob("*.whl"))
    sdists = list(dist.glob("*.tar.gz"))
    assert len(wheels) == len(sdists) == 1, "Expected exactly one wheel and one sdist"
    assert {p.name for p in dist.iterdir()} == {wheels[0].name, sdists[0].name}, "Unexpected artifact"
    license_text = (project / "LICENSE").read_bytes()
    assert b"MIT License" in license_text
    assert b"Copyright (c) 2026 panyw5" in license_text

    with ZipFile(wheels[0]) as archive:
        names = archive.namelist()
        for sub in SUBPACKAGES:
            assert f"pyw/{sub}/__init__.py" in names, f"Missing subpackage: {sub}"
        license_name = next(n for n in names if n.endswith("/LICENSE"))
        assert archive.read(license_name) == license_text
        metadata = archive.read(next(n for n in names if n.endswith("/METADATA"))).decode()
        assert "License-File: LICENSE" in metadata
        assert f"Version: {version}\n" in metadata
        metadata_name = next(line[6:] for line in metadata.splitlines() if line.startswith("Name: "))
        normalize = lambda name: re.sub(r"[-_.]+", "-", name).lower()
        assert normalize(metadata_name) == normalize(distribution_name), "Wrong distribution name"
        assert archive.read("pyw/core/root_system.py") == (project / "pyw/core/root_system.py").read_bytes()
        check_names(names)
    with tarfile.open(sdists[0]) as archive:
        names = archive.getnames()
        license_name = next(n for n in names if n.endswith("/LICENSE"))
        member = archive.extractfile(license_name)
        assert member is not None and member.read() == license_text
        check_names(names)

    import pyw
    from pyw.core.root_system import AffineRootSystem

    installed = Path(pyw.__file__).resolve()
    assert not installed.is_relative_to(project.resolve()), "Imported source checkout, not installed wheel"
    assert pyw.__version__ == version
    print(f"Installed package: {installed}", flush=True)
    blocks = re.findall(r"```python\n(.*?)```", (project / "README.md").read_text(), re.S)
    assert len(blocks) == 2, "Review verifier when changing README examples"
    for i, block in enumerate(blocks, 1):
        exec(compile(block, f"README-example-{i}", "exec"), {})
        print(f"README example {i}: PASS", flush=True)
    for ct, count in [(["A", 2], 6), (["B", 2], 8), (["G", 2], 12)]:
        assert len(AffineRootSystem(ct).coroots) == count
    try:
        AffineRootSystem(["A", 2, 1]).coroots
    except NotImplementedError:
        pass
    else:
        raise AssertionError("Affine coroot enumeration must fail explicitly")
    print("Artifacts, license, versions and installed coroot interface: PASS", flush=True)

    if run_tests:
        import pytest

        # Copy only the reviewed tests, not the pyw package; test imports must
        # resolve to the installed wheel rather than to the source checkout.
        with tempfile.TemporaryDirectory(prefix="pyw-release-tests-") as directory:
            target = Path(directory)
            for name in (*TEST_FILES, "conftest.py"):
                shutil.copyfile(project / "pyw/tests" / name, target / name)
            result = pytest.main([
                "-c", "/dev/null", "-o", "addopts=",
                "-o", "markers=sage: requires SageMath", "-p", "no:cacheprovider",
                "--import-mode=importlib", "-q", "--tb=short", str(target),
            ])
            if result != 0:
                raise SystemExit(int(result))


def check_names(names: list[str]) -> None:
    forbidden = {"tests", "__pycache__", ".git", ".env", "config.json"}
    for name in names:
        assert not forbidden.intersection(Path(name).parts), f"Unexpected file: {name}"
        assert not name.endswith((".ipynb", ".pyc", ".DS_Store")), f"Unexpected file: {name}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--dist", type=Path, required=True)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    verify(args.project.resolve(), args.dist.resolve(), args.run_tests)
