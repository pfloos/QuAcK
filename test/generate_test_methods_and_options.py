import re

enabled = set()

# Load methods that should be enabled
with open("methods.enable") as f:
    for line in f:
        name = line.strip()
        if name:
            enabled.add(name)


def transform_line(line):
    out = []
    for token in line.split():
        if token in ("T", "F"):
            # This is a flag, skip: replaced later
            continue
        out.append(token)
    return out


with open("../input/methods.default") as f_in, open("methods.test", "w") as f_out:
    for line in f_in:
        if line.strip().startswith("#"):
            # Header line, just copy
            f_out.write(line)
            # Parse method names on this header line
            current_methods = transform_line(line[1:])
            continue

        if not line.strip():
            f_out.write(line)
            continue

        # This is the line with flags: rebuild it
        tokens = line.split()
        flags = []
        for method in current_methods:
            flags.append("T" if method in enabled else "F")

        f_out.write("  " + "   ".join(flags) + "\n")


# --- Load overrides --------------------------------------------------------

overrides = {}
current_block = None

with open("options.override") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        # A header line like "# HF:" or "# spin:"
        if line.startswith("#"):
            block = line[1:].strip().rstrip(":")
            current_block = block
            overrides.setdefault(current_block, {})
            continue

        # A line like "triplet = T"
        if "=" in line and current_block is not None:
            key, value = [x.strip() for x in line.split("=", 1)]
            overrides[current_block][key] = value


# --- Process default options ------------------------------------------------

def parse_header(line):
    """Extract block name and option names."""
    m = re.match(r"#\s*(.*?):\s*(.*)", line)
    if not m:
        return None, None
    block_name = m.group(1).strip()
    options = m.group(2).split()
    return block_name, options


with open("../input/options.default") as fin, open("options.test", "w") as fout:
    current_block = None
    current_keys = None

    for line in fin:
        if line.strip().startswith("#"):
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
