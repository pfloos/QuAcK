import numpy as np


def generate_cap_file(filename, size, width, seed=42):
    """
    Generates a file with random CAP integral values in the format:
    mu nu wx wy wz
    """
    np.random.seed(seed)  # For reproducibility
    with open(filename, 'w') as f:
        for mu in range(1, size + 1):
            for nu in range(mu, size + 1):  # Only upper triangle to avoid duplicate entries
                # Generate three random values
                wx, wy, wz = np.random.rand(3)*width
                f.write(f"{mu} {nu} {wx:.6f} {wy:.6f} {wz:.6f}\n")


def read_and_construct_matrix(filename, size):
    """
    Reads the file and constructs the symmetric matrix W.
    """
    W = np.zeros((size, size))
    with open(filename, 'r') as f:
        for line in f:
            mu, nu, wx, wy, wz = line.split()
            mu, nu = int(mu) - 1, int(nu) - 1  # Convert to zero-based index
            value = float(wx) + float(wy) + float(wz)
            W[mu, nu] = value
            W[nu, mu] = value  # Enforce symmetry
    return W


if __name__ == "__main__":
    with open("../int/nBas.dat", 'r') as f:
        size = int(f.readline().strip())
    print("nBas: ", size)
    width = 0
    generate_cap_file("../int/CAP.dat", size, width)
    W = read_and_construct_matrix("../int/CAP.dat", size)
    print("W matrix:")
    print(W)
