
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
from utils import print_col


current_date = datetime.now()

quack_root = os.getenv('QUACK_ROOT')

# User Name
user_name = os.getlogin()

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



print(f"The current date and time is {current_date.strftime('%Y-%m-%d %H:%M:%S')}")
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

args = parser.parse_args()

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

    def __init__(self, mol, multip, basis, geom, methd):
        self.mol = mol
        self.multip = multip
        self.basis = basis
        self.geom = geom
        self.methd = methd

    def prep_inp(self):

        # geometry
        generate_xyz(self.geom, filename="{}.xyz".format(self.mol))

        # input files
        for inp in ["methods", "options"]:
            inp_file = "{}.{}".format(inp, self.methd.upper())
            if os.path.exists("inp/{}".format(inp_file)):
                shutil.copy("inp/{}".format(inp_file), "../mol/{}".format(inp_file))
            else:
                print_col("File 'inp/{}' does not exist.".format(inp_file), "red")
                sys.exit(1)

    def run(file_out, mol, bas, multip):

        os.chdir('..')
        print(f" :$ cd ..")

        for file_in in ["methods", "options"]:
            command = ['cp', 'tests/{}.RHF'.format(file_in), 'input/{}'.format(file_in)]
            print(f" :$ {' '.join(command)}")
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode != 0:
                print("Error moving file: {}".format(result.stderr))

        command = [
            'python{}'.format(PYTHON_VERSION), 'PyDuck.py',
            '-x', '{}'.format(mol), 
            '-b', '{}'.format(bas),
            '-m', '{}'.format(multip)
        ]
        print(f" :$ {' '.join(command)}")
        with open(file_out, 'w') as fobj:
            result = subprocess.run(command, stdout=fobj, stderr=subprocess.PIPE, text=True)
        if result.stderr:
            print("Error output:", result.stderr)

        os.chdir('tests')
        print(f" :$ cd tests")


# ---


def main():

    work_path = Path('{}/tests/work'.format(quack_root))
    if not work_path.exists():
        work_path.mkdir(parents=True, exist_ok=True)
        print(f"Directory '{work_path}' created.\n")

    for mol in molecules:

        mol_name = mol.name
        mol_mult = mol.multiplicity
        mol_geom = mol.geometry
        mol_data = mol.properties

        print_col("  Molecule: {} (2S+1 = {})".format(mol_name, mol_mult), "blue")

        for mol_prop_name, mol_prop_data in mol_data.items():

            print_col("    Testing {}".format(mol_prop_name), "cyan")

            methd = mol_prop_name[len('properties_'):]

            if(len(mol_prop_data) == 0):
                print_col("    {} is empty. Skipping...".format(mol_prop_name), "cyan")
                print()
                continue

            for basis_name, basis_data in mol_prop_data.items():
                print_col("      Basis set = {}".format(basis_name), "yellow")

                if(len(basis_data) == 0):
                    print_col("      {} is empty. Skipping...".format(basis_name), "yellow")
                    print()
                    continue

                work_methd = Path('{}/{}'.format(work_path, methd))
                if not work_methd.exists():
                    work_methd.mkdir(parents=True, exist_ok=True)
                    #print(f"Directory '{work_methd}' created.\n")
        
                New_Quack_Job = Quack_Job(mol_name, mol_mult, basis_name, mol_geom, methd)
                New_Quack_Job.prep_inp()

#                for name, val in basis_data.items():
#                    print(f"      name = {name}")
#                    print(f"      val = {val}")

                print()
            print()
        print()

    quit()

        
#                # create input files
#                class_methd.gen_input()
#        
#                file_out = "{}/{}/{}_{}_{}.out".format(work_path, prop, mol_name, mol_mult, bas)
#
#                print(" testing {} for {}@{} (2S+1 = {})".format(prop, mol_name, bas, mol_mult))
#                print(" file_out: {}".format(file_out))
#
#                class_methd.run_job(file_out, mol_name, bas, mol_mult)


        


db_name = '{}.db'.format(bench)

molecules = get_molecules_from_db(db_name)

main()

