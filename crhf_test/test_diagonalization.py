import numpy as np
import re


def read_complex_matrix(filename, N):
    with open(filename, "r") as file:
        lines = file.readlines()

    # Flattened list of complex numbers
    complex_numbers = []

    # Regular expression to match complex numbers in the format (a,b)
    complex_pattern = re.compile(
        r"\(\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?),\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\)")

    for line in lines:
        matches = complex_pattern.findall(line)
        for real, imag in matches:
            complex_numbers.append(complex(float(real), float(imag)))

    # Reshape into NxN matrix
    return np.array(complex_numbers, dtype=np.complex128).reshape(N, N)


# Read the matrix stored from the first core guess and compare values to the fortran code --> wokrs
N = 10  # Set matrix size
matrix = read_complex_matrix("matrix_test", N)
print("Matrix")
print(matrix)
print("Element 5,1")
print(matrix[4, 0])
eigenval, eigenvec = np.linalg.eig(matrix)
print("Eigenvalues")
print(eigenval)
print("Eigenvectors")
print(eigenvec)
