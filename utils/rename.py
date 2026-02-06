#!/usr/bin/env python3
import sys
from pathlib import Path

def main():
    if len(sys.argv) != 2:
        print("1 argument required!!")
        sys.exit(1)

    suffix = sys.argv[1]

    files_to_rename = ["Ov.dat", "Kin.dat", "Nuc.dat", "ERI.dat", "F12.dat", "Erf.dat"]

    for fname in files_to_rename:
        src = Path(fname)
        if src.exists():
            dst = Path(f"{src.stem}.{suffix}{src.suffix}")
            src.rename(dst)
            print(f"Renamed {src} → {dst}")
        else:
            print(f"Warning: {src} does not exist, skipping")

if __name__ == "__main__":
    main()
