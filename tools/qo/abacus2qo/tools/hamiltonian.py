import numpy as np
from abacus2qo.tools.basic_functions import make_complex

def recover_from_upper_triangle(matrix: list) -> np.ndarray:
    """recover full matrix from upper-right triangle

    Args:
        matrix (list): matrix of upper-right triangle

    Returns:
        np.ndarray: full matrix
    """
    matrix = np.array(matrix)
    matrix = matrix + matrix.conj().T
    for i in range(len(matrix)):
        matrix[i][i] = matrix[i][i] / 2
    return matrix

def parse(nkpts: int, path = "./") -> tuple[list, list, int]:
    """parse hamiltonian and overlap matrix H(k) and S(k)

    Args:
        nkpts (int): number of kpoints
    """
    hamiltonian = []
    overlap = []
    num_bands = 0
    if path[-1] != "/":
        path += "/"
    for ik in range(nkpts):
        h_fname = path + f"data-{ik}-H"
        s_fname = path + f"data-{ik}-S"
        
        hamiltonian_k = []
        overlap_k = []
        # first read Hamiltonian H(k)
        with open(h_fname, "r") as f:
            lines = f.readlines()
        for il, line in enumerate(lines):
            if len(line.strip()) == 0:
                continue
            hamiltonian_k_b = []
            words = line.strip().split()
            if il == 0:
                num_bands = int(words[0])
                # delete the first element
                words.pop(0)
            for icol in range(il):
                hamiltonian_k_b.append(complex(0, 0))
            for word in words:
                hamiltonian_k_b.append(make_complex(word))

            hamiltonian_k.append(hamiltonian_k_b)

        hamiltonian_k = recover_from_upper_triangle(hamiltonian_k)
        hamiltonian.append(hamiltonian_k)

        # then read overlap matrix S(k)
        with open(s_fname, "r") as f:
            lines = f.readlines()
        for il, line in enumerate(lines):
            if len(line.strip()) == 0:
                continue
            overlap_k_b = []
            words = line.strip().split()
            if il == 0:
                num_bands = int(words[0])
                # delete the first element
                words.pop(0)
            for icol in range(il):
                overlap_k_b.append(complex(0, 0))
            for word in words:
                overlap_k_b.append(make_complex(word))

            overlap_k.append(overlap_k_b)

        overlap_k = recover_from_upper_triangle(overlap_k)
        overlap.append(overlap_k)

    #return recover_from_upper_triangle(np.array(hamiltonian)), recover_from_upper_triangle(np.array(overlap)), num_bands
    return hamiltonian, overlap, num_bands

if __name__ == "__main__":
    pass