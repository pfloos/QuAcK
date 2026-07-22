#!/usr/bin/env python3

from pathlib import Path
import subprocess


README = Path("../README.md")
PYDUCK = Path("../PyDuck.py")


def get_pyduck_help():
    """Run PyDuck.py -h and return its output."""
    result = subprocess.run(
        ["python", str(PYDUCK), "-h"],
        capture_output=True,
        text=True,
        check=True
    )

    return result.stdout


def get_options():
    """Run PyDuck.py -h and return its output."""
    result = subprocess.run(
        ["cat", "../input/options.default"],
        capture_output=True,
        text=True,
        check=True
    )

    return result.stdout


def get_methods():
    """Run PyDuck.py -h and return its output."""
    result = subprocess.run(
        ["cat", "../input/methods.default"],
        capture_output=True,
        text=True,
        check=True
    )

    return result.stdout


def update_section(text, name, content):
    """Replace content between BEGIN/END markers."""

    begin = f"<!-- BEGIN {name} -->"
    end = f"<!-- END {name} -->"

    if begin not in text or end not in text:
        raise RuntimeError(f"Missing markers for {name}")

    before = text.split(begin)[0]
    after = text.split(end)[1]

    return (
        before
        + begin
        + '\n'
        + content.rstrip()
        + '\n'
        + end
        + after
    )


def main():

    readme = README.read_text()

    help_text = get_pyduck_help()
    options_text = get_options()
    methods_text = get_methods()
    readme = update_section(
        readme,
        "PYDUCK_HELP",
        help_text
    )
    readme = update_section(
        readme,
        "options",
        options_text
    )
    readme = update_section(
        readme,
        "methods",
        methods_text
    )

    README.write_text(readme)


if __name__ == "__main__":
    main()
