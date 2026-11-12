"""Run the H2 examples from the QuAcK manual."""

import argparse
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = {
    "rhf": {"RHF"},
    "g0w0": {"RHF", "G0W0"},
    "ftda": {"RHF", "FTDA"},
    "water-ftda": {"RHF", "FTDA"},
}


def write_methods(path, enabled):
    lines = (ROOT / "input/methods.default").read_text().splitlines()
    seen = set()
    for index, line in enumerate(lines[:-1]):
        if line.startswith("#"):
            flags = line[1:].split()
            lines[index + 1] = " ".join("T" if flag in enabled else "F" for flag in flags)
            seen.update(flags)
    if enabled - seen:
        raise ValueError(f"Example flags missing from template: {enabled - seen}")
    path.write_text("\n".join(lines) + "\n")


def check_result(example, output):
    if "Total wall time for QuAcK" not in output:
        raise ValueError("QuAcK did not reach the final timing line")
    expressions = {
        "rhf": r"RHF\s+energy\s*=\s*(\S+)",
        "g0w0": r"G0W0@RHF HOMO\s+energy\s*=\s*(\S+)",
        "ftda": r"1p weight Z, residual \(Ha\)\s*\n\s*1\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",
    }
    if example == "water-ftda":
        _, found, table = output.partition("1p weight Z, residual (Ha)")
        rows = [line.split() for line in table.splitlines()[1:5]]
        if not found or len(rows) != 4 or any(
            len(row) != 5 or row[0] != str(index)
            or not all(math.isfinite(float(value)) for value in row[1:])
            for index, row in enumerate(rows, 1)
        ):
            raise ValueError("Expected four finite FTDA roots for water")
        return
    match = re.search(expressions[example], output)
    if match is None:
        raise ValueError(f"No result found for {example}")
    if not all(math.isfinite(float(field)) for field in match.groups()):
        raise ValueError(f"Nonfinite result for {example}")


def run(example, directory, installation):
    (directory / "input").mkdir()
    (directory / "mol").mkdir()
    write_methods(directory / "input/methods", EXAMPLES[example])
    options = (ROOT / "input/options.default").read_text()
    molecule = "water" if example == "water-ftda" else "H2"
    if example == "water-ftda":
        options, count = re.subn(
            r"(^# Faddeev:[^\n]*\n[^\n]*?)[0-9]+[ \t]*$", r"\g<1>4", options, flags=re.M
        )
        if count != 1:
            raise ValueError("Expected one standalone Faddeev options record")
    (directory / "input/options").write_text(options)
    geometry = f"mol/{molecule}.xyz"
    shutil.copyfile(ROOT / geometry, directory / geometry)
    env = dict(os.environ, QUACK_ROOT=str(installation))
    result = subprocess.run(
        [
            sys.executable, str(ROOT / "PyDuck.py"), "--working_dir", str(directory),
            "-x", molecule, "-b", str(ROOT / "basis/6-31g"), "--no_cap",
        ],
        cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    print(result.stdout, end="")
    result.check_returncode()
    check_result(example, result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("example", choices=EXAMPLES)
    parser.add_argument("--work-dir", type=Path, help="New directory to retain the inputs")
    args = parser.parse_args()
    installation = Path(os.environ.get("QUACK_ROOT", ROOT)).resolve()
    if not (installation / "bin/QuAcK").is_file():
        parser.error(f"Build QuAcK first. Missing {installation / 'bin/QuAcK'}")
    if args.work_dir:
        directory = args.work_dir.resolve()
        directory.mkdir()
        run(args.example, directory, installation)
    else:
        with TemporaryDirectory(prefix="quack-docs-") as temporary:
            run(args.example, Path(temporary), installation)


if __name__ == "__main__":
    main()
