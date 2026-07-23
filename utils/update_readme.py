from pathlib import Path
import subprocess
import os


QUACK_ROOT = Path(
    os.environ.get(
        "QUACK_ROOT",
        Path(__file__).resolve().parent.parent
    )
)

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
    """Return the content of the options.default file."""
    options_text = (
        QUACK_ROOT / "input/options.default"
    ).read_text()
    print(options_text)
    return f"""```
QuAcK 💩 % cat input/options
{options_text}
```"""


def get_methods():
    """Return the content of the methods.default file."""
    methods_text = (
        QUACK_ROOT / "input/methods.default"
    ).read_text()

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
        + "\n"
        + content.rstrip()
        + "\n"
        + end
        + after
    )


def update_readme(readme, sections):
    """
    Update all generated README sections.
    """

    for section in sections:
        section_text = globals()[f"get_{section}"]()
        readme = update_section(
            readme,
            section,
            section_text
        )

    return readme


def main():

    readme = README.read_text()

    sections = [
        "pyduck_help",
        "options",
        "methods",
    ]

    readme = update_readme(
        readme,
        sections
    )

    README.write_text(readme)


if __name__ == "__main__":
    main()
