
import argparse

from molecule import save_molecules_to_json, load_molecules_from_json
from molecule import create_database, add_molecule_to_db, remove_database

from feather_bench import FeatherBench


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


db_name = '{}.db'.format(bench)


# Save molecules to JSON
#save_molecules_to_json(FeatherBench, 'FeatherBench.json')

# Load molecules from JSON
#loaded_molecules = load_molecules_from_json('FeatherBench.json')
#print(loaded_molecules)

#remove_database(db_name)

create_database(db_name)
for molecule in FeatherBench:
    add_molecule_to_db(db_name, molecule)



