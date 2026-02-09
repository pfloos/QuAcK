#!/usr/bin/env python3
import sys
from pathlib import Path

def usage():
    print("Usage: create_subroutine.py <SubroutineName>")
    sys.exit(1)

def main():
    if len(sys.argv) != 2:
        usage()

    name = sys.argv[1].strip()
    fname = f"{name}.f90"
    path = Path(fname)

    if path.exists():
        print(f"Error: {fname} already exists.")
        sys.exit(1)

    print(f"Creating {fname}...")

    content = f"""! Automatically generated subroutine {name}
subroutine {name}()
    implicit none

    ! TODO: implement {name}

end subroutine {name}
"""

    with open(path, "w") as f:
        f.write(content)

    print("Done!")

if __name__ == "__main__":
    main()
