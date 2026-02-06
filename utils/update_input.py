from pathlib import Path
import shutil

def copy_file(src, dst):
    print(f"copying {src} to {dst}")
    shutil.copy(src, dst)

def main():
    copy_file(Path("input/methods.default"), Path("input/methods"))
    copy_file(Path("input/options.default"), Path("input/options"))
    print("done!")

if __name__ == "__main__":
    main()
