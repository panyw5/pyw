# Release checklist

pyw is an Alpha research library requiring SageMath. Do not upload a release
until the checks below have completed and the maintainer has approved it.

## Known release-preparation limitations

- The maintainer approved the MIT license with copyright holder
  `2026 panyw5`. The full text is in [LICENSE](../LICENSE); verify that both
  distribution formats include it whenever rebuilding release candidates.
- The `coroots` property returns all positive and negative coroots for finite
  types. For affine (or other non-finite) types it now explicitly raises
  `NotImplementedError` rather than silently returning only a generating basis.
  Use `simple_coroots` when a generating basis is required. This is an intentional
  correction to the previous misleading behavior.
- Full-suite coverage and the complete character test suite are not yet
  verified by the minimum release checks. A targeted passing test selection
  must not be presented as a passing full suite.
- Adjacent `test_principal.py` and `test_nondegenerate.py` checks produced
  10 passes and 15 failures in SageMath 10.4 / Python 3.12.4. These include
  missing public attributes, Python Fraction/Sage coercion, Weyl action domain
  mismatches, and invalid finite-type test indices. They require separate
  implementation-versus-test triage before those features are advertised as
  verified. The packaging and core fractional fixes do not resolve them.
  The maintainer has deferred this work: these incomplete experimental APIs
  are outside the minimum release scope, remain in the package, and must not
  be advertised as validated functionality. Their failures remain recorded
  rather than being skipped or removed from the test suite.

## Before building

- Confirm the version in `pyproject.toml` matches `pyw.__version__`.
- Confirm the license text, copyright holder, and rights to all distributed code.
- Review the README's supported features, examples, environment requirements,
  and project links. Do not advertise unimplemented features as available.
- Run the test suite in the intended SageMath environment. Record failures,
  expected failures, skipped tests, and coverage, not just the exit status.
- Confirm that no credentials, private research data, caches, notebooks, or
  development-only files are included in the distribution.

## Build and validate

Use a clean checkout or a temporary copy so build artifacts do not obscure
uncommitted changes. Build tools do not need SageMath to generate the artifacts;
runtime validation does.

```bash
python -m pip install build twine
python -m build
python -m twine check --strict dist/*
```

Inspect both the wheel and sdist. Required functional subpackages must be
present. The sdist must be sufficient to rebuild a working wheel.

## Installed-package checks

Create a temporary environment using SageMath's Python, with access to its
installed modules. Install the built wheel without resolving SageMath again:

```bash
sage -python -m venv --system-site-packages /tmp/pyw-release-venv
/tmp/pyw-release-venv/bin/python -m pip install --no-deps dist/pyw_algebra-<version>-py3-none-any.whl
```

Change to a directory outside the source checkout and check `import pyw`, its
version and module path, and both README examples. Use a unique temporary
environment directory rather than overwriting an existing environment.

Also test installation from the sdist. Verify imports come from the temporary
installation, not an editable checkout. Validate ordinary pip dependency
resolution separately if promising installation without preinstalled SageMath.

## Publication

Only after explicit approval, upload to TestPyPI, validate the installed
artifacts there, then publish to PyPI. These are separate services. A 404 for a
project name does not guarantee that the name is available for registration.
Do not commit tokens to the repository. Prefer trusted publishing or use
credentials managed outside the source tree.

PyPI distribution filenames cannot be reused for replacement uploads. Treat
each published version as immutable and use a new version for corrections.

### GitHub OIDC trusted publishing

The default publishing path is [publish.yml](../.github/workflows/publish.yml).
It is manually dispatched on `main`, defaults to TestPyPI, and has three jobs:

1. **build**: isolated PEP 517 build (wheel rebuilt from sdist), strict Twine
   metadata checks, and upload of a single immutable run artifact.
2. **validate**: SageMath 10.4 container installs the wheel and executes
   [verify_release.py](../scripts/verify_release.py). This checks source-external
   imports, README examples, versions, licenses, package contents, coroot behavior,
   and the six reviewed test files against the installed package. It deliberately
   does not claim that the incomplete principal/nondegenerate APIs pass tests.
   The SageMath 10.4 image is run as `linux/amd64`, matching the GitHub Ubuntu
   runner. Local Apple Silicon reproduction requires amd64 emulation; this image
   does not provide an arm64 manifest.
3. **publish**: downloads exactly the artifact validated above and uploads using
   a short-lived OIDC credential. This job has `id-token: write`, does not check
   out source, and does not run the package. It also requests publish attestations.

No API token or `TWINE_PASSWORD` secret is needed. Twine remains a local/build
metadata checker, not the default uploader.

#### Required site-side setup (not performed by editing repository files)

Create **both** GitHub Environments under repository Settings -> Environments:

| Environment | Purpose |
| --- | --- |
| `testpypi` | TestPyPI uploads |
| `pypi` | Production uploads |

Restrict deployment branches to `main`. Configure required reviewers and prevent
self-review where your GitHub plan/repository supports it; protect `main` and
review changes to the workflow and verifier. Merely referencing an Environment
in YAML does **not** enable approval protection. If approval rules are unavailable,
do not claim the workflow has an enforced human approval gate.

Configure a GitHub Trusted Publisher independently on TestPyPI and PyPI. If `pyw-algebra`
does not yet exist on a service, use that service's **pending publisher** setup;
if it already exists and you own it, use its project Publishing settings.

| Field | Exact value |
| --- | --- |
| PyPI project name (pending publisher) | `pyw-algebra` |
| GitHub owner | `panyw5` |
| Repository name | `pyw` |
| Workflow filename | `publish.yml` (not the full path) |
| Environment name on TestPyPI | `testpypi` |
| Environment name on PyPI | `pypi` |

Do not authorize a different workflow or omit the Environment binding. A
publisher configured on one service does not authorize uploads to the other.

#### Trigger and verify

Commit and push the reviewed release files and workflow only with maintainer
approval. In GitHub Actions, select **Publish Python package (OIDC)** ->
**Run workflow**, choose branch `main`, and leave destination `testpypi` for
the first run. Review build and Sage validation results before approving the
Environment deployment. Running or approving this workflow is an upload action.

After TestPyPI installation checks pass, explicitly approve a separate run with
destination `pypi`. Each run rebuilds and revalidates its selected commit; publish
the same reviewed commit to both services. A failed or partial upload is not
silently skipped: inspect the uploaded files and version before retrying.

The initial configuration uses major/release Action tags and the SageMath 10.4
image tag. For stronger supply-chain controls, pin reviewed Actions to full
commit SHAs and the container to a verified digest, and maintain those pins.
GitHub-hosted container execution and OIDC issuance must be verified on the
first workflow run; local Sage tests do not verify GitHub's environment.

Once uploaded, install into a fresh temporary SageMath-compatible environment:

```bash
python -m pip install --no-deps --index-url https://test.pypi.org/simple/ pyw-algebra==0.1.0
```

Here `python` must be the temporary environment's SageMath-compatible Python.
Do not combine TestPyPI and PyPI indexes for this check: SageMath is already
provided, and `--no-deps` avoids resolving test-index dependencies. Repeat the
source-external import, README, and coroot checks on the downloaded package.
