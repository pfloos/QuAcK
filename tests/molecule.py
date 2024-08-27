
import json
import sqlite3

class Molecule:
    def __init__(self, name, multiplicity, geometry, energies):
        self.name = name
        self.multiplicity = multiplicity
        self.geometry = geometry  # List of tuples (atom, x, y, z)
        self.energies = energies  # Dictionary of dictionaries: {method: {basis: energy}}

    def get_energy(self, method, basis):
        """Retrieve energy for a specific method and basis set."""
        return self.energies.get(method, {}).get(basis, None)

    def to_dict(self):
        return {
            "name": self.name,
            "multiplicity": self.multiplicity,
            "geometry": self.geometry,
            "energies": self.energies,
        }

    @staticmethod
    def from_dict(data):
        return Molecule(
            name=data["name"],
            multiplicity=data["multiplicity"],
            geometry=data["geometry"],
            energies=data["energies"]
        )

def save_molecules_to_json(molecules, filename):
    with open(filename, 'w') as f:
        json_data = [molecule.to_dict() for molecule in molecules]
        json.dump(json_data, f, indent=4)

def load_molecules_from_json(filename):
    with open(filename, 'r') as f:
        json_data = json.load(f)
        return [Molecule.from_dict(data) for data in json_data]


def create_database(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS molecules
                      (name TEXT, multiplicity INTEGER, geometry TEXT, energies TEXT)''')
    conn.commit()
    conn.close()

def add_molecule_to_db(db_name, molecule):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    geometry_str = json.dumps(molecule.geometry)
    energies_str = json.dumps(molecule.energies)
    cursor.execute("INSERT INTO molecules VALUES (?, ?, ?, ?)",
                   (molecule.name, molecule.multiplicity, geometry_str, energies_str))
    conn.commit()
    conn.close()

def get_molecules_from_db(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute("SELECT name, multiplicity, geometry, energies FROM molecules")
    rows = cursor.fetchall()
    molecules = []
    for row in rows:
        name, multiplicity, geometry_str, energies_str = row
        geometry = json.loads(geometry_str)
        energies = json.loads(energies_str)  # energies is a dictionary of dictionaries
        molecules.append(Molecule(name, multiplicity, geometry, energies))
    conn.close()
    return molecules

def write_geometry_to_xyz(molecule, filename):
    with open(filename, 'w') as f:
        # First line: number of atoms
        f.write(f"{len(molecule.geometry)}\n")
        # Second line: empty comment line
        f.write("\n")
        # Remaining lines: atom positions
        for atom, x, y, z in molecule.geometry:
            f.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")


