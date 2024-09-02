
import sys


def read_quantities_from_file(filename):
    quantities = {}

    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 2):
            # Remove any leading or trailing whitespace/newline characters
            quantity_name = lines[i].strip()
            quantity_value = float(lines[i+1].strip())
            quantities[quantity_name] = quantity_value

    return quantities

def print_quantities(quantities):
    for key, value in quantities.items():
        print(f'"{key}": {value},')

filename = sys.argv[1]

quantities = read_quantities_from_file(filename)
print_quantities(quantities)

