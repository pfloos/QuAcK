
import time
import threading
import sys
import os
import shutil
from pathlib import Path
import subprocess
import platform
from datetime import datetime
import argparse

from molecule import get_molecules_from_db
from molecule import generate_xyz
from utils import print_col, stdout_col


current_date = datetime.now()

quack_root = os.getenv('QUACK_ROOT')

# User Name
# Fallback to 'default_user' if USER is not set
user_name = os.getenv('USER', 'default_user')

# Operating System
os_name = platform.system()
os_release = platform.release()
os_version = platform.version()

# CPU Information
machine = platform.machine()
processor = platform.processor()

# System Architecture
architecture = platform.architecture()[0]

# Python Version
python_version_full = platform.python_version_tuple()
PYTHON_VERSION = "{}.{}".format(python_version_full[0], python_version_full[1])


print(f"The current date and time is {
      current_date.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"User Name: {user_name}")
print(f"Operating System: {os_name} {os_release} ({os_version})")
print(f"CPU: {processor} ({machine})")
print(f"System Architecture: {architecture}")
print(f"QUACK_ROOT: {quack_root}")
print(f"Python version: {python_version_full}\n\n")


parser = argparse.ArgumentParser(description="Benchmark Data Sets")

parser.add_argument(
    '-s', '--set_type',
    choices=['light', 'medium', 'heavy'],
    default='light',
    help="Specify the type of data set: light (default), medium, or heavy."
)

thresh_default = 1e-7
parser.add_argument(
    '-t', '--thresh',
    type=float,
    default=thresh_default,
    help='Threshold for acceptable difference, default = {}'.format(
        thresh_default)
)


args = parser.parse_args()

THRESH = args.thresh

if args.set_type == 'light':
    bench = 'FeatherBench'
    bench_title = "\n\nSelected Light Benchmark: {}\n\n".format(bench)
elif args.set_type == 'medium':
    bench = 'BalanceBench'
    bench_title = "\n\nSelected Medium Benchmark: {}\n\n".format(bench)
elif args.set_type == 'heavy':
    bench = 'TitanBench'
    bench_title = "\n\nSelected Heavy Benchmark: {}\n\n".format(bench)
else:
    bench_title = "\n\nSelected Light Benchmark: {}\n\n".format(bench)

print(bench_title.center(150, '-'))
print("\n\n")

# ---


class Quack_Job:

    def __init__(self, mol, multip, basis, geom, methd, workdir):
        self.mol = mol
        self.multip = multip
        self.basis = basis
        self.geom = geom
        self.methd = methd
        self.workdir = workdir

    def prep_inp(self):

        # geometry
        if not os.path.exists("{}/mol".format(self.workdir)):
            os.makedirs("{}/mol".format(self.workdir))
        generate_xyz(
            self.geom, filename="{}/mol/{}.xyz".format(self.workdir, self.mol))

        # input files
        for inp in ["methods", "options"]:
            inp_file = "{}.{}".format(inp, self.methd.upper())
            if os.path.exists("inp/{}".format(inp_file)):
                shutil.copy("{}/tests/inp/{}".format(quack_root, inp_file),
                            "{}/input/{}".format(self.workdir, inp))
            else:
                print_col(
                    "File 'inp/{}' does not exist.".format(inp_file), "red")
                sys.exit(1)

    def run_ci(self):

        try:

            os.chdir('..')

            command = [
                'python{}'.format(PYTHON_VERSION), 'PyDuck.py',
                '--working_dir', '{}'.format(self.workdir),
                '-x', '{}'.format(self.mol),
                '-b', '{}'.format(self.basis),
                '-m', '{}'.format(self.multip),
                '--no_cap'
            ]

            file_out = "{}/{}/{}_{}_{}.out".format(
                self.workdir, self.methd, self.mol, self.multip, self.basis)
            with open(file_out, 'w') as fobj:
                result = subprocess.run(
                    command, stdout=fobj, stderr=subprocess.PIPE, text=True)
            if result.stderr:
                print("Error output:", result.stderr)
                print("QuAcK output of failed test:")
                print(result.stdout)
                sys.exit(1)

            os.chdir('tests')

        except Exception as e:

            print_col(f"An error occurred: {str(e)}", "red")

    def run(self):

        def display_spinner():
            spinner = ['|', '/', '-', '\\']
            idx = 0
            while not done_event.is_set():
                stdout_col(f'\r    Testing {self.methd} ({self.basis}) {
                           spinner[idx]}', "cyan")
                sys.stdout.flush()
                idx = (idx + 1) % len(spinner)
                time.sleep(0.05)
            stdout_col(f'\r    Testing {self.methd} ({
                       self.basis})    \n\n', "cyan")

        done_event = threading.Event()
        spinner_thread = threading.Thread(target=display_spinner)
        spinner_thread.start()

        try:

            os.chdir('..')
            # print_col(f"      Starting QuAck..", "magenta")
            # print_col(f"      $ cd ..", "magenta")

            command = [
                'python{}'.format(PYTHON_VERSION), 'PyDuck.py',
                '--working_dir', '{}'.format(self.workdir),
                '-x', '{}'.format(self.mol),
                '-b', '{}'.format(self.basis),
                '-m', '{}'.format(self.multip),
                '--no_cap'
            ]
            # print_col(f"      $ {' '.join(command)}", "magenta")

            file_out = "{}/{}/{}_{}_{}.out".format(
                self.workdir, self.methd, self.mol, self.multip, self.basis)
            with open(file_out, 'w') as fobj:
                result = subprocess.run(
                    command, stdout=fobj, stderr=subprocess.PIPE, text=True)
            if result.stderr:
                print("Error output:", result.stderr)
                sys.exit(1)

            os.chdir('tests')
            # print_col(f"      $ cd tests", "magenta")

        except Exception as e:

            print_col(f"An error occurred: {str(e)}", "red")

        finally:

            done_event.set()
            spinner_thread.join()

    def check_data(self, data_ref, test_failed_):
        filepath = '{}/test/Rtest.dat'.format(self.workdir)
        data_new = {}
        try:
            # read data_new
            with open(filepath, 'r') as f:
                lines = f.readlines()
                for i in range(0, len(lines) - 1, 2):
                    key = lines[i].strip()
                    value = lines[i + 1].strip()
                    data_new[key] = float(value)  # Convert value to float

            # Compare with data_ref
            for key in data_ref:
                if key not in data_new:
                    print_col(f"        😐 {key} missing ⚠️ ", "yellow")
                    test_failed_ = True
                else:
                    diff = abs(data_new[key] - data_ref[key]
                               ) / (1e-15 + abs(data_ref[key]))
                    if (diff <= THRESH):
                        print_col(f"        🙂 {key}", "green")
                    else:
                        print_col(f"        ☹️  {key}: ❌ {
                                  data_ref[key]} ≠ {data_new[key]}", "red")
                        test_failed_ = True
        except FileNotFoundError:
            print_col(f"Error: The file '{filepath}' does not exist.", "red")
            sys.exit(1)
        except Exception as e:
            print_col(f"An error occurred: {str(e)}", "red")
            sys.exit(1)


# ---


def main():

    work_path = Path('{}/tests/work'.format(quack_root))
    if not work_path.exists():
        work_path.mkdir(parents=True, exist_ok=True)
        print(f"Directory '{work_path}' created.\n")

    # for I/O
    if not os.path.exists("{}/test".format(work_path)):
        os.makedirs("{}/test".format(work_path))
    if not os.path.exists("{}/input".format(work_path)):
        os.makedirs("{}/input".format(work_path))

    is_ci = os.getenv('CI')
    test_failed = False
    for mol in molecules:

        mol_name = mol.name
        mol_mult = mol.multiplicity
        mol_geom = mol.geometry
        mol_data = mol.properties

        print_col("  Molecule: {} (2S+1 = {})".format(mol_name, mol_mult), "blue")

        for mol_prop_name, mol_prop_data in mol_data.items():

            methd = mol_prop_name[len('properties_'):]

            if (len(mol_prop_data) == 0):
                continue

            for basis_name, basis_data in mol_prop_data.items():

                if (len(basis_data) == 0):
                    continue

                work_methd = Path('{}/{}'.format(work_path, methd))
                if not work_methd.exists():
                    work_methd.mkdir(parents=True, exist_ok=True)

                New_Quack_Job = Quack_Job(
                    mol_name, mol_mult, basis_name, mol_geom, methd, work_path)
                New_Quack_Job.prep_inp()
                if is_ci:
                    New_Quack_Job.run_ci()
                else:
                    New_Quack_Job.run()

                test_failed_ = False
                New_Quack_Job.check_data(basis_data, test_failed_)
                if (test_failed_):
                    test_failed = True

                print()
            print()
        print()

    if test_failed:
        sys.exit(1)

    sys.exit(0)


db_name = '{}.db'.format(bench)

molecules = get_molecules_from_db(db_name)

main()
