import numpy as np

# Function to parse a complex matrix from a file
def parse_complex_matrix(file_path):
    matrix = []
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line into elements
            elements = line.strip().split()
            row = []
            for element in elements:
                # Remove parentheses and split into real and imaginary parts
                real, imag = element.strip('()').split(',')
                # Convert to complex number
                row.append(complex(float(real), float(imag)))
            matrix.append(row)
    return np.array(matrix)

# Function to write a complex matrix to a file
def write_complex_matrix(matrix, file_path):
    with open(file_path, 'w') as file:
        for row in matrix:
            # Format each element as (x,y)
            formatted_row = ' '.join(f"({z.real},{z.imag})" for z in row)
            file.write(formatted_row + '\n')

# File paths
einsum_file = 'einsum.log'
gemm_file = 'gemm.log'
output_file = 'difference.log'

# Parse matrices from files
matrix1 = parse_complex_matrix(einsum_file)
matrix2 = parse_complex_matrix(gemm_file)

# Check if matrices have the same shape
if matrix1.shape != matrix2.shape:
    print("Error: Matrices must have the same shape.")
else:
    # Compute the difference
    difference_matrix = matrix1 - matrix2

    # Write the difference matrix to a new file
    write_complex_matrix(difference_matrix, output_file)

    # Compute the modulus of the difference matrix
    modulus_matrix = np.abs(difference_matrix)

    # Find the maximum modulus value
    max_modulus = np.max(modulus_matrix)

    # Output the result
    print(f"The maximum modulus of the difference is: {max_modulus}")