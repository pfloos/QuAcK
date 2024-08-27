
import sqlite3

from molecule import Molecule
from molecule import save_molecules_to_json, load_molecules_from_json
from molecule import create_database, add_molecule_to_db



molecules = [
    Molecule(
        name="H2O",
        multiplicity=1,
        geometry=[
            {"element": "O", "x": 0.000000, "y": 0.000000, "z": 0.117790},
            {"element": "H", "x": 0.000000, "y": 0.755453, "z": -0.471161},
            {"element": "H", "x": 0.000000, "y": -0.755453, "z": -0.471161}
        ],
        energies={
            "RHF": {
                "cc-pvdz": -76.0267058009, 
                "cc-pvtz": -76.0570239304, 
                "cc-pvqz": -76.0646816616
            },
        }
    ),
]


# Save molecules to JSON
save_molecules_to_json(molecules, 'molecules.json')

# Load molecules from JSON
loaded_molecules = load_molecules_from_json('molecules.json')
print(loaded_molecules)

# Create a database and add molecules
db_name = 'molecules.db'
create_database(db_name)
for molecule in molecules:
    add_molecule_to_db(db_name, molecule)

