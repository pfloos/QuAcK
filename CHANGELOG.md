# Changelog

## 1.0 (2026-11-12)

QuAcK’s first stable source release brings together the work since the 0.9 prerelease.

- The method inventory covers restricted, unrestricted, and generalized Hartree–Fock reference paths where implemented, with post-Hartree–Fock methods including coupled cluster, RPA, Green’s-function methods, and Bethe–Salpeter response. See the [method availability table](https://pfloos.github.io/QuAcK/generated/methods.html) for scope by reference.
- The canonical RHF driver offers [single-shot GW–Parquet](docs/gw-parquet.md) and separate [FTDA and FRPA workflows](docs/faddeev.md). For calculations using the older Parquet/Faddeev inputs, use the bundled templates and follow the migration guidance in the standalone guide.
- The PyDuck front end uses PySCF for molecular integrals. The [quickstart](docs/quickstart.md) includes runnable examples for H₂ RHF, single-shot GW, and FTDA, and for water ionization with single-shot GW.
- The [manual](docs/manual.md) documents the input flags, their implementation coverage, and the source routines. The repository includes a light benchmark suite and checks for the documented examples.
- Calculation output reports the software version and, for a Git checkout, the source commit. [Citation metadata](CITATION.cff) records the release version. The [README](README.md#-contributors) lists the contributors.

[Compare 0.9 and v1.0](https://github.com/pfloos/QuAcK/compare/0.9...v1.0)

## 0.9 (2020-04-09)

[Original GitHub prerelease](https://github.com/pfloos/QuAcK/releases/tag/0.9)
