#ifndef CSR_READER_H
#define CSR_READER_H

#include <fstream>

#include "file_reader.h"
#include "sparse_matrix.h"

namespace ModuleIO
{

/**
 * @brief Class to read CSR file
 * @details This class is used to read CSR file
 * the default format is like:
 * ```
 *  --- Ionic Step 1 ---
 *  # print density matrix in real space DM(R)
 *  26 # number of localized basis
 *  13 # number of Bravais lattice vector R
 *
 *  unitcell_information
 *
 *  #--------------------------------------------------------#
 *  #                      CSR format                        #
 *  # Outer loop is the number of Bravais lattice vectors.   #
 *  # First line is Bravais lattice vector index Rx, Ry, Rz, #
 *  # followed by the number of non-zero elements.           #
 *  # Next are three blocks of data.                         #
 *  #--------------------------------------------------------#
 * 
 *  -1 0 0 507
 *  # CSR values
 *  6.73361941e-04 -3.97537783e-05 7.92408228e-04 ...
 *  # CSR column indices
 *  0 1 2 ...
 *  # CSR row pointers
 *  0 26 52 ...
 *  ...
 *
 * ```
 * It will store the R coordinates and sparse matrices as two vectors.
 * One can use getter functions to get the R coordinates and sparse matrices,
 * and related info including step, matrix dimension, number of R.
 */
template <typename T>
class csrFileReader : public FileReader
{
  public:
    // Constructor
    csrFileReader(const std::string& filename);

    // read all matrices of all R coordinates
    void parseFile();

    // get number of R
    int getNumberOfR() const;

    // get sparse matrix of a specific R coordinate
    SparseMatrix<T> getMatrix(int Rx, int Ry, int Rz) const;

    // get matrix by using index
    SparseMatrix<T> getMatrix(int index) const;

    // get R coordinate using index
    std::vector<int> getRCoordinate(int index) const;

    // get step
    int getStep() const;

    // get matrix dimension
    int getMatrixDimension() const;

  private:
    std::vector<std::vector<int>> RCoordinates;
    std::vector<SparseMatrix<T>> sparse_matrices;
    int step;
    int matrixDimension;
    int numberOfR;
};

} // namespace ModuleIO

#endif // READ_CSR_H
