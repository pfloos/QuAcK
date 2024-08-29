
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

    def run(self, work_path):

        def display_spinner():
            spinner = ['|', '/', '-', '\\']
            idx = 0
            while not done_event.is_set():
                stdout_col(f'\r    Testing {self.methd} ({self.basis}) {spinner[idx]}', "yellow")
                sys.stdout.flush()
                idx = (idx + 1) % len(spinner)
                time.sleep(0.1)
            stdout_col(f'\r    Testing {self.methd} ({self.basis})    ', "yellow")

        done_event = threading.Event()
        spinner_thread = threading.Thread(target=display_spinner)
        spinner_thread.start()

        try:
    
            os.chdir('..')
            #print_col(f"      Starting QuAck..", "magenta")
            #print_col(f"      $ cd ..", "magenta")
    
            command = [
                'python{}'.format(PYTHON_VERSION), 'PyDuck.py',
                '-x', '{}'.format(self.mol), 
                '-b', '{}'.format(self.basis),
                '-m', '{}'.format(self.multip)
            ]
            #print_col(f"      $ {' '.join(command)}", "magenta")
    
            file_out = "{}/{}/{}_{}_{}.out".format(work_path, self.methd, self.mol, self.multip, self.basis)
            with open(file_out, 'w') as fobj:
                result = subprocess.run(command, stdout=fobj, stderr=subprocess.PIPE, text=True)
            if result.stderr:
                print("Error output:", result.stderr)
    
            os.chdir('tests')
            #print_col(f"      $ cd tests", "magenta")
    
        except Exception as e:

            print_col(f"An error occurred: {str(e)}", "red")

        finally:

            done_event.set()
            spinner_thread.join()
    
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

            #print_col("    Testing {}".format(mol_prop_name), "cyan")

            methd = mol_prop_name[len('properties_'):]

            if(len(mol_prop_data) == 0):
                #print_col("    {} is empty. Skipping...".format(mol_prop_name), "cyan")
                #print()
                continue

            for basis_name, basis_data in mol_prop_data.items():
                #print_col("      Basis set: {}".format(basis_name), "yellow")

                if(len(basis_data) == 0):
                    #print_col("      {} is empty. Skipping...".format(basis_name), "yellow")
                    #print()
                    continue

                work_methd = Path('{}/{}'.format(work_path, methd))
                if not work_methd.exists():
                    work_methd.mkdir(parents=True, exist_ok=True)
                    #print(f"Directory '{work_methd}' created.\n")
        
                New_Quack_Job = Quack_Job(mol_name, mol_mult, basis_name, mol_geom, methd)
                New_Quack_Job.prep_inp()
                New_Quack_Job.run(work_path)

#                for name, val in basis_data.items():
#                    print(f"      name = {name}")
#                    print(f"      val = {val}")

                print()
            print()
        print()

    quit()

        
db_name = '{}.db'.format(bench)

molecules = get_molecules_from_db(db_name)

main()

