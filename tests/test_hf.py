
import os
from pathlib import Path
import subprocess
import platform
from datetime import datetime

from molecule import get_molecules_from_db


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

# ---

mp2 = "# MP2 MP3\n  F   F\n"
cc = "# CCD pCCD DCD CCSD CCSD(T)\n  F   F    F   F    F\n"
rcc = "# drCCD rCCD crCCD lCCD\n  F     F    F     F\n"
ci = "# CIS CIS(D) CID CISD FCI\n  F   F      F   F    F\n"
rpa = "# phRPA phRPAx crRPA ppRPA\n  F     F      F     F\n"
gf = "# G0F2 evGF2 qsGF2 ufGF2 G0F3 evGF3\n  F    F     F     F    F    F\n"
gw = "# G0W0 evGW qsGW SRG-qsGW ufG0W0 ufGW\n  F    F    F    F        F      F\n"
gtpp = "# G0T0pp evGTpp qsGTpp ufG0T0pp\n  F      F      F      F\n"
gteh = "# G0T0eh evGTeh qsGTeh\n  F      F      F\n"
tests = "# Rtest Utest Gtest\n  F     F     F\n"

# ---

hf_opt = "# HF: maxSCF thresh  DIIS guess mix shift stab search\n      256    0.00001 5    1     0.0 0.0   F    F\n"
mp_opt = "# MP: reg\n      F\n"
cc_opt = "# CC: maxSCF thresh   DIIS\n      64     0.00001  5\n"
tda_opt = "# spin: TDA singlet triplet\n        F   T       T\n"
gf_opt = "# GF: maxSCF thresh  DIIS lin eta renorm reg\n      256    0.00001 5    F   0.0 0      F\n"
gw_opt = "# GW: maxSCF thresh  DIIS lin eta TDA_W reg\n      256    0.00001 5    F   0.0 F     F\n"
gt_opt = "# GT: maxSCF thresh  DIIS lin eta TDA_T reg\n      256    0.00001 5    F   0.0 F     F\n"
acfdt_opt = "# ACFDT: AC Kx  XBS\n         F  F   T\n"
bse_opt = "# BSE: phBSE phBSE2 ppBSE dBSE dTDA\n        F     F      F     F    T\n"
list_opt = [hf_opt, mp_opt, cc_opt, tda_opt, gf_opt, gw_opt, gt_opt, acfdt_opt, bse_opt]

# ---

class class_RHF:

    def gen_input():

        f = open("methods", "w")
        f.write("# RHF UHF GHF ROHF\n")
        f.write("  T   F   F   F\n")
        f.write("{}{}{}{}{}{}{}{}{}{}".format(mp2, cc, rcc, ci, rpa, gf, gw, gtpp, gteh, tests))
        f.close()

        f = open("options", "w")
        for opt in list_opt:
            f.write("{}".format(opt))
        f.close()

    def run_job(file_out, mol, bas, multip):

        os.chdir('..')
        print(f" :$ cd ..")

        for file_in in ["methods", "options"]:
            command = ['cp', 'tests/{}'.format(file_in), 'input/{}'.format(file_in)]
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

class class_UHF:
    def gen_input():
        f = open("methods", "w")
        f.write("# RHF UHF GHF ROHF\n")
        f.write("  F   T   F   F\n")
        f.write("{}{}{}{}{}{}{}{}{}{}".format(mp2, cc, rcc, ci, rpa, gf, gw, gtpp, gteh, tests))
        f.close()

# ---

class class_GHF:
    def gen_input():
        f = open("methods", "w")
        f.write("# RHF UHF GHF ROHF\n")
        f.write("  F   F   T   F\n")
        f.write("{}{}{}{}{}{}{}{}{}{}".format(mp2, cc, rcc, ci, rpa, gf, gw, gtpp, gteh, tests))
        f.close()

# ---

class class_ROHF:
    def gen_input():
        f = open("methods", "w")
        f.write("# RHF UHF GHF ROHF\n")
        f.write("  F   F   F   T\n")
        f.write("{}{}{}{}{}{}{}{}{}{}".format(mp2, cc, rcc, ci, rpa, gf, gw, gtpp, gteh, tests))
        f.close()

# ---

class_map = {
    "RHF": class_RHF,
    "UHF": class_UHF,
    "GHF": class_GHF,
    "ROHF": class_ROHF,
}

def main():

    work_path = Path('{}/tests/work'.format(quack_root))
    if not work_path.exists():
        work_path.mkdir(parents=True, exist_ok=True)
        print(f"Directory '{work_path}' created.\n")

    for mol in molecules:

        mol_name = mol.name
        mol_mult = mol.multiplicity


        for methd in list_methd:

            if methd not in mol.energies:
                print(f"Method {methd} does not exist for {mol_name}.")
                continue

            for bas, _ in mol.energies[methd].items():

                work_methd = Path('{}/{}'.format(work_path, methd))
                if not work_methd.exists():
                    work_methd.mkdir(parents=True, exist_ok=True)
                    print(f"Directory '{work_methd}' created.\n")
        
                class_methd = class_map.get(methd)
        
                # create input files
                class_methd.gen_input()
        
                file_out = "{}/{}/{}_{}_{}.out".format(work_path, methd, mol_name, mol_mult, bas)

                print(" testing {} for {}@{} (2S+1 = {})".format(methd, mol_name, bas, mol_mult))
                print(" file_out: {}".format(file_out))

                class_methd.run_job(file_out, mol_name, bas, mol_mult)

                print("\n")
            print("\n\n")

        print(" --- --- --- ---")
        print("\n\n\n")
        

db_name = 'molecules.db'
molecules = get_molecules_from_db(db_name)

list_methd = ["RHF", "UHF", "GHF", "ROHF"]

main()

