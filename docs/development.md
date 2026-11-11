# Validation and documentation

## Find the implementation

PyDuck writes molecular and integral files. `src/QuAcK/QuAcK.f90` reads the calculation inputs and dispatches to the R, U, G, or B driver. The drivers then call the method-specific directories in `src/`. Use the [Fortran source index](generated/routines.md) to locate individual files and procedure names.

## Run the checks

The regular CI builds QuAcK and runs the light benchmark suite. Its [workflow](https://github.com/pfloos/QuAcK/blob/master/.github/workflows/ci.yml) shows the database and fixture setup. Method-specific scripts in `tests/` check GW–Parquet and Faddeev calculations. Read each script's arguments and fixtures before running it.

The three [documentation examples](quickstart.md) are smoke tests. They check that the selected drivers produce finite results and that QuAcK finishes. Benchmark tests are the place for numerical reference comparisons.

## Build the manual locally

From the repository root, with the PyDuck dependencies installed:

```bash
python -m pip install -r docs/requirements.txt
python utils/build_docs_reference.py
python -m mkdocs build --strict
python -m mkdocs serve
```

The reference generator reads the CLI help, shipped templates, curated method inventory, and tracked Fortran sources. It checks that every method flag has a row and that all inventory source links resolve in the checkout. Generated pages carry the Git revision. Keep the generated files out of commits.

When adding a flag, update the input reader and template, then record the active R/U/G paths and limitations in `docs/methods.tsv`. For a change to an existing flag, review its row and the relevant method guide. Where possible, record or extend a small test that exercises the normal entry point and checks an observable.

## Publish

The documentation workflow builds on pull requests and pushes to `master`. To publish from the repository, choose **GitHub Actions** under **Settings → Pages → Build and deployment**. A build from `master` then deploys the manual from its own commit. The build also provides an HTML artifact for review.
