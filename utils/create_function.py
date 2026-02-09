#!/usr/bin/env python3
import sys
from pathlib import Path

def usage():
    print("Usage: create_function.py <FunctionName> [ReturnType]")
    print("Example: create_function.py compute_energy real(dp)")
    sys.exit(1)

def main():
    if len(sys.argv) < 2:
        usage()

    name = sys.argv[1].strip()
    return_type = sys.argv[2] if len(sys.argv) > 2 else "real(dp)"

    fname = f"{name}.f90"
    path = Path(fname)

    if path.exists():
        print(f"Error: {fname} already exists.")
        sys.exit(1)

    print(f"Creating {fname}...")

    content = f"""! Automatically generated function {name}
{return_type} function {name}()
    implicit none

    ! TODO: implement {name}

    {name} = 0.0
end function {name}
"""

    with open(path, "w") as f:
        f.write(content)

    print("Done!")

if __name__ == "__main__":
    main()
