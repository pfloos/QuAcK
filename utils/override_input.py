import os
import re


def parse_header(line):
    """Extract block name and option names."""
    m = re.match(r"#\s*(.*?):\s*(.*)", line)
    if not m:
        return None, None
    block_name = m.group(1).strip()
    options = m.group(2).split()
    return block_name, options


def override_options(overrides, options_file):
    """
    Overrides specified options in the option file.

    Parameters
    ----------
    overrides   : Dictionary containing the option as key and the value as value
                  with which its value should be replaced.
    options_file : Filepath of the options file. 

    """

    # --- Process options ------------------------------------------------

    with open(options_file) as fin, open('options.tmp', 'w') as fout:
        current_block = None
        current_keys = None

        for line in fin:
            if line.strip().startswith('#'):
                # Write header line unchanged
                fout.write(line)

                # Parse new block
                block, keys = parse_header(line)
                current_block = block
                current_keys = keys
                continue

            # Value line
            if current_block and current_keys and line.strip():

                values = line.split()

                # Apply overrides if any
                if current_block in overrides:
                    for i, key in enumerate(current_keys):
                        if key in overrides[current_block]:
                            values[i] = overrides[current_block][key]

                # Reconstruct line
                fout.write("      " + "   ".join(values) + "\n")
                continue

            # Blank lines / formatting untouched
            fout.write(line)

    os.replace('options.tmp', options_file)
