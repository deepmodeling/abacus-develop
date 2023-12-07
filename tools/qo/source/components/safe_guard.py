import numpy as np

class toQO_SafeGuard:
    """I am not that soap. Use to check hermiticity, orthonormality, rank, linear dependence of matrices
    """
    def __init__(self) -> None:
        pass

    def initialize(self) -> None:
        """initialize the class
        """
        pass

    def check_hermiticity(self, mat: np.ndarray) -> float:
        """check if a matrix is Hermite

        Args:
            mat (np.ndarray): matrix to be checked
        """
        error = np.linalg.norm(mat - mat.conj().T)
        print("Error of matrix being Hermite is: %14.10e"%error)
        return error

    def check_orthonormality(self, mat: np.ndarray) -> float:
        """check if a matrix is orthonormal

        Args:
            mat (np.ndarray): matrix to be checked
        """
        error = np.linalg.norm(mat - np.eye(mat.shape[-1]))
        print("Error of matrix being orthonormal is: %14.10e"%error)
        return error

    def check_identity(self, mat: np.ndarray) -> float:
        """check if a matrix is identity

        Args:
            mat (np.ndarray): matrix to be checked
        """
        error = np.linalg.norm(mat - np.eye(mat.shape[-1]))
        print("Error of matrix being identity is: %14.10e"%error)
        return error

    def check_diagonal(self, mat: np.ndarray) -> float:
        """check if a matrix is diagonal

        Args:
            mat (np.ndarray): matrix to be checked
        """
        error = np.linalg.norm(mat - np.diag(np.diag(mat)))
        print("Error of matrix being diagonal is: %14.10e"%error)
        return error

    def check_rank(self, mat: np.ndarray) -> int:
        """check rank of a matrix

        Args:
            mat (np.ndarray): matrix to be checked
        """
        rank = np.linalg.matrix_rank(mat)
        print("Rank of matrix is: ", rank)
        return rank

    def check_linear_dependence(self, mat: np.ndarray) -> int:
        """check if a matrix is linearly dependent

        Args:
            mat (np.ndarray): matrix to be checked
        """
        rank = np.linalg.matrix_rank(mat)
        if rank < mat.shape[-1]:
            print("Matrix is linearly dependent!")
            exit(1)
        else:
            return rank

    def check_eigenstates(self, eigenvec: np.ndarray, ovlp: np.ndarray) -> float:
        """check if eigenstates are orthonormal

        Args:
            eigenvec (np.ndarray): eigenstates
            ovlp (np.ndarray): overlap matrix
        """
        error = np.linalg.norm(eigenvec.conj().T @ ovlp @ eigenvec - np.eye(eigenvec.shape[-1]))
        print("Error of eigenstates being orthonormal is: %14.10e"%error)

        return error