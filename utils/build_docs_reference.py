"""Build the input and method reference from the current QuAcK checkout."""

import csv
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "docs" / "generated"
REPOSITORY = "https://github.com/pfloos/QuAcK"


def git_revision():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def source_link(path, revision):
    target = ROOT / path
    if not target.exists():
        raise ValueError(f"Missing documentation source path: {path}")
    kind = "tree" if target.is_dir() else "blob"
    return f"[`{path}`]({REPOSITORY}/{kind}/{revision}/{path})"


def code_block(value):
    return f"```text\n{value.rstrip()}\n```\n"


def build_inputs(revision):
    env = dict(os.environ, QUACK_ROOT=str(ROOT), COLUMNS="88")
    help_text = subprocess.check_output(
        [sys.executable, "PyDuck.py", "--help"], cwd=ROOT, env=env, text=True
    )
    help_text = help_text[help_text.index("usage:"):]
    parts = [
        "# Input reference",
        f"Source revision `{revision}`",
        "## PyDuck command line",
        code_block(help_text),
        "Set `QUACK_ROOT` before running PyDuck. The charge parser takes an integer, "
        "so use `--charge=-1` for an anion.",
    ]
    for name in ("methods", "options"):
        template = f"input/{name}.default"
        reader = f"src/QuAcK/read_{name}.f90"
        parts.extend([
            f"## `input/{name}`",
            f"Shipped template {source_link(template, revision)}",
            f"Fortran reader {source_link(reader, revision)}",
            code_block((ROOT / template).read_text()),
        ])
    return "\n\n".join(parts)


def build_methods(revision):
    lines = [
        "# Method availability",
        "",
        f"Source revision `{revision}`",
        "",
        "**R** = RHF, **U** = UHF, **G** = GHF in the standard drivers. "
        "**yes** marks an active source path and **—** marks an unavailable path. "
        "**B path** refers to the separate Bogoliubov driver. "
        "Source coverage and numerical validation are tracked separately.",
        "",
        "| Method or flag | R | U | G | Scope | Source |",
        "|---|---|---|---|---|---|",
    ]
    expected = {
        flag for line in (ROOT / "input/methods.default").read_text().splitlines()
        if line.startswith("#") for flag in line[1:].split()
    }
    covered = set()
    with (ROOT / "docs/methods.tsv").open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            flags = row["Method or input flag"].split(", ")
            duplicates = covered.intersection(flags)
            if duplicates:
                raise ValueError(f"Duplicate method flags: {sorted(duplicates)}")
            covered.update(flags)
            labels = ", ".join(f"`{flag}`" for flag in flags)
            source = source_link(row["Source"], revision)
            lines.append(
                f"| {labels} | {row['R']} | {row['U']} | {row['G']} | "
                f"{row['Scope']} | {source} |"
            )
    if covered != expected:
        raise ValueError(f"Missing flags: {sorted(expected - covered)}. "
                         f"Obsolete flags: {sorted(covered - expected)}.")
    lines.extend([
        "",
        "### Other references",
        "",
        "ROHF uses the restricted driver. Its post-HF coverage must be checked "
        "method by method. Complex RHF and UHF have dedicated response paths. "
        "The Bogoliubov driver is selected with HFB or BRPA and has its own "
        "branches for MP2, response, scGF2 and scGW.",
        "",
        "### Compatibility",
        "",
        "FTDA and FRPA use the [standalone workflow](../faddeev.md). "
        "For pre-separation inputs see the "
        "[historical implementation notes](../TDA_Faddeev_Parquet.md).",
    ])
    return "\n".join(lines)


def build_routines(revision):
    lines = ["# Fortran source index", "", f"Source revision `{revision}`", ""]
    pattern = re.compile(r"^\s*(?:\w+\s+)*(subroutine|function)\s+(\w+)\b", re.I)
    tracked = subprocess.check_output(
        ["git", "ls-files", "-z", "src"], cwd=ROOT, text=True
    )
    sources = {}
    for name in tracked.split("\0"):
        path = ROOT / name
        if path.suffix.lower() == ".f90":
            sources.setdefault(path.parent, []).append(path)
    for folder, files in sorted(sources.items()):
        lines.extend([f"## {folder.relative_to(ROOT / 'src')}", ""])
        for path in sorted(files):
            names = []
            for line in path.read_text(errors="replace").splitlines():
                match = pattern.match(line)
                if match and not line.lstrip().lower().startswith("end"):
                    names.append(match[2])
            link = source_link(path.relative_to(ROOT).as_posix(), revision)
            routines = ", ".join(f"`{name}`" for name in names)
            lines.append(f"- {link}" + (f" — {routines}" if routines else ""))
        lines.append("")
    return "\n".join(lines)


def main():
    revision = git_revision()
    pages = {
        "inputs.md": build_inputs(revision),
        "methods.md": build_methods(revision),
        "routines.md": build_routines(revision),
    }
    GENERATED.mkdir(parents=True, exist_ok=True)
    for filename, contents in pages.items():
        (GENERATED / filename).write_text(contents.rstrip() + "\n")


if __name__ == "__main__":
    main()
