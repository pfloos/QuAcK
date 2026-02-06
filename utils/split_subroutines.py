import re
from pathlib import Path
import shutil

SOURCE_DIR = Path(".")

subroutine_start_re = re.compile(
    r'^\s*(?:recursive\s+|pure\s+|elemental\s+|impure\s+)*subroutine\s+(\w+)',
    re.IGNORECASE
)

subroutine_end_re = re.compile(
    r'^\s*end\s+subroutine\b',
    re.IGNORECASE
)


def extract_subroutines(lines):
    """Return list of (name, full_text) for each subroutine."""
    subroutines = []
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]
        start_match = subroutine_start_re.match(line)

        if start_match:
            name = start_match.group(1)
            block_lines = [line]
            i += 1
            depth = 1

            while i < n and depth > 0:
                current_line = lines[i]

                if subroutine_start_re.match(current_line):
                    depth += 1
                if subroutine_end_re.match(current_line):
                    depth -= 1

                block_lines.append(current_line)
                i += 1

            subroutines.append((name, "".join(block_lines)))
        else:
            i += 1

    return subroutines


def main():
    for filepath in SOURCE_DIR.glob("*.f90"):
        print(f"\nProcessing {filepath.name}")

        with open(filepath, "r", encoding="utf-8") as f:
            lines = f.readlines()

        subs = extract_subroutines(lines)

        if len(subs) <= 1:
            print("  0 or 1 subroutine — leaving untouched.")
            continue

        print(f"  Found {len(subs)} subroutines → splitting (build-safe mode)")

        # Backup original
        backup_path = filepath.with_suffix(filepath.suffix + ".old")
        shutil.move(filepath, backup_path)
        print(f"  Original backed up as {backup_path.name}")

        # First subroutine stays under original filename
        first_name, first_text = subs[0]
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(first_text)
        print(f"  -> {filepath.name} now contains only: {first_name}")

        # Remaining subroutines get their own files
        for name, text in subs[1:]:
            out_file = SOURCE_DIR / f"{name}.f90"
            print(f"  -> Writing {out_file.name}")
            with open(out_file, "w", encoding="utf-8") as f:
                f.write(text)


if __name__ == "__main__":
    main()
