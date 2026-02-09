import re
from pathlib import Path

ROOT_DIR = Path(".")

subroutine_start_re = re.compile(
    r'^\s*(?:recursive\s+|pure\s+|elemental\s+|impure\s+)*subroutine\s+(\w+)',
    re.IGNORECASE
)


def count_subroutines(filepath):
    count = 0
    names = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            match = subroutine_start_re.match(line)
            if match:
                count += 1
                names.append(match.group(1))

    return count, names


def main():
    print("Scanning for files with multiple subroutines...\n")

    found_any = False

    for filepath in ROOT_DIR.rglob("*.f90"):
        count, names = count_subroutines(filepath)

        if count > 1:
            found_any = True
            rel = filepath.relative_to(ROOT_DIR)
            print(f"{rel}  →  {count} subroutines")
            print("   ", ", ".join(names))
            print()

    if not found_any:
        print("No files contain more than one subroutine.")


if __name__ == "__main__":
    main()
