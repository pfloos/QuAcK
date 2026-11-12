# Install and run

## Build

Install GNU Fortran, a C compiler, BLAS/LAPACK, HDF5 development files, Make, Ninja, Git, and Python. On Ubuntu, the project CI also installs `wget`, `m4`, and `pkg-config`. See the [CI build](https://github.com/pfloos/QuAcK/blob/master/.github/workflows/ci.yml) for the Ubuntu package list. On macOS, install the equivalent dependencies and follow the compiler and library settings in [the build generator](https://github.com/pfloos/QuAcK/blob/master/src/make_ninja.py).

```bash
git clone https://github.com/pfloos/QuAcK.git
cd QuAcK
export QUACK_ROOT="$PWD"
python3 -m venv .venv
source .venv/bin/activate
python -m pip install pyscf
make -C src
```

For the stable release, clone with `git clone --branch v1.0 https://github.com/pfloos/QuAcK.git`. Keep the checkout path free of spaces for the build.

Use `./configure -i ninja` from the repository root if Ninja needs to be bootstrapped. The build produces `bin/QuAcK`. CAP calculations also need PyOpenCAP. The examples below disable CAP integrals.

## Run H₂

The [example runner](examples/run.py) prepares a fresh calculation directory from the shipped input templates, uses the bundled H₂ geometry and 6-31G basis, and writes the output to standard output.

```bash
python docs/examples/run.py rhf
python docs/examples/run.py g0w0
```

To keep the input, integrals, and output together, pass a new directory and redirect the log.

```bash
python docs/examples/run.py g0w0 --work-dir h2-gw > h2-gw.log
```

To set up a calculation by hand, copy `input/methods.default` to `input/methods` and `input/options.default` to `input/options` in the calculation directory. Enable the reference and target method by changing the corresponding `F` values to `T`. Keep the records in order. Place the geometry under `mol/` and pass its name without `.xyz` to `-x`. Use an absolute basis-file path or a name under the calculation's `basis/` directory.

```bash
python "$QUACK_ROOT/PyDuck.py" --working_dir "$PWD" -b 6-31g -x H2 --no_cap
```

## Read the output

The bundled H₂ geometry is at 1.235 Å. With the shipped 6-31G basis and options,
the local examples give the following values.

| Example | Observable | Value |
|---|---|---:|
| RHF | Total energy | −1.0487553440 Ha |
| Single-shot GW | HOMO quasiparticle energy | −13.499288 eV |


- RHF prints the total energy and HF orbital energies. Read the labels for units.
- Single-shot GW prints quasiparticle energies and renormalization factors. For a removal pole referenced to vacuum, the ionization energy is the negative of the pole energy.

Check the SCF and any outer iterations, look for a completed QuAcK timing line, and check that the requested method printed results. PyDuck launches the Fortran program as a child process, so its exit status alone is insufficient to establish that the calculation finished.

## Ionize water

Run `python docs/examples/run.py water-gw` from the repository root.
The example uses `mol/water.xyz` and the 6-31G basis. For occupied orbitals,
the vertical ionization energy is the negative of the reported `e_GW`.
The figure shows the four valence quasiparticles and their weights `Z`.

| RHF orbital | IP (eV) | Quasiparticle weight |
|---|---:|---:|
| 5 | 12.05 | 0.9569 |
| 4 | 14.00 | 0.9562 |
| 3 | 18.70 | 0.9618 |
| 2 | 34.87 | 0.6737 |
