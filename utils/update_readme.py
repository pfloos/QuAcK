#!/usr/bin/env python3

from pathlib import Path
import subprocess
import os

QUACK_ROOT = Path(os.environ.get("QUACK_ROOT", ".."))

README = QUACK_ROOT / "README.md"
PYDUCK = QUACK_ROOT / "PyDuck.py"


def get_pyduck_help():
    """Run PyDuck.py -h and return its output."""
    help_text = subprocess.run(
        ["python", str(PYDUCK), "-h"],
        capture_output=True,
        text=True,
        check=True
    ).stdout
    return f"""```
~ 💩 % cd $QUACK_ROOT
QuAcK 💩 % python PyDuck.py -h
{help_text.rstrip()}
```"""


def get_options():
    """Return the content of the options file output."""
    options_text = subprocess.run(
        ["cat", QUACK_ROOT / "input/options.default"],
        capture_output=True,
        text=True,
        check=True
    ).stdout
    return f"""```
QuAcK 💩 % cat input/options 
{options_text.rstrip()}
```"""


def get_methods():
    """Return the content of the methods file."""
    methods_text = subprocess.run(
        ["cat", QUACK_ROOT / "input/methods.default"],
        capture_output=True,
        text=True,
        check=True
    ).stdout
    return f"""```
QuAcK 💩 % cat input/methods 
{methods_text.rstrip()}
```"""


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


def update_readme(readme, sections):
    """
    Updates all sections in the readme text by the text defined by the
    get_%section function

    Parameters
    ----------
    readme : Text of the Readme
    sections :  Sections which are supposed to be updated

    Returns
    -------

    """
    for section in sections:
        section_text = globals()[f'get_{section}']()
        readme = update_section(readme, section, section_text)
    return readme


def main():

    readme = README.read_text()
    sections = ['pyduck_help', 'options', 'methods']
    update_readme(readme, sections)
    README.write_text(readme)


if __name__ == "__main__":
    main()
