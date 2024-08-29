
import sqlite3

from molecule import save_molecules_to_json, load_molecules_from_json
from molecule import create_database, add_molecule_to_db

from swift_set import swiftset



# Save molecules to JSON
save_molecules_to_json(swiftset, 'swiftset.json')

# Load molecules from JSON
loaded_molecules = load_molecules_from_json('swiftset.json')
print(loaded_molecules)

# Create a database and add molecules
db_name = 'swiftset.db'
create_database(db_name)
for molecule in swiftset:
    add_molecule_to_db(db_name, molecule)



