#include "csr_reader.h"

#include "module_base/tool_quit.h"

namespace ModuleIO
{

// constructor
template <typename T>
csrFileReader<T>::csrFileReader(const std::string& filename) : FileReader(filename)
{
    parseFile();
}

// function to parse file
template <typename T>
void csrFileReader<T>::parseFile()
{
    // Check if file is open
    if (!isOpen())
    {
        ModuleBase::WARNING_QUIT("csrFileReader::parseFile", "File is not open");
    }

    std::string tmp_string;

    // Read the step
    readLine();
    ss >> tmp_string >> step;

    // Read the matrix dimension
    readLine();
    ss >> tmp_string >> tmp_string >> tmp_string >> tmp_string >> matrixDimension;

    // Read the number of R
    readLine();
    ss >> tmp_string >> tmp_string >> tmp_string >> tmp_string >> numberOfR;

    // Read the matrices
    for (int i = 0; i < numberOfR; i++)
    {
        std::vector<int> RCoord(3);
        int nonZero;

        readLine();
        ss >> RCoord[0] >> RCoord[1] >> RCoord[2] >> nonZero;
        RCoordinates.push_back(RCoord);
        if (nonZero > matrixDimension * matrixDimension)
        {
            ModuleBase::WARNING_QUIT("csrFileReader::parseFile",
                                     "Number of non-zero elements is larger than matrix dimension");
        }
        SparseMatrix<T> matrix(matrixDimension, matrixDimension);
        if (nonZero > 0)
        {

            std::vector<T> csr_values(nonZero);
            std::vector<int> csr_col_ind(nonZero);
            std::vector<int> csr_row_ptr(matrixDimension + 1);

            // read CSR values
            readLine();
            for (int i = 0; i < nonZero; i++)
            {
                ss >> csr_values[i];
            }
            // read column indices
            readLine();
            for (int i = 0; i < nonZero; i++)
            {
                ss >> csr_col_ind[i];
            }
            // read row pointers
            readLine();
            for (int i = 0; i < matrixDimension + 1; i++)
            {
                ss >> csr_row_ptr[i];
            }

            matrix.readCSR(csr_values, csr_col_ind, csr_row_ptr);
        }
        sparse_matrices.push_back(matrix);
    }
}

// function to get R coordinate
template <typename T>
std::vector<int> csrFileReader<T>::getRCoordinate(int index) const
{
    if (index < 0 || index >= RCoordinates.size())
    {
        ModuleBase::WARNING_QUIT("csrFileReader::getRCoordinate", "Index out of range");
    }
    return RCoordinates[index];
}

// function to get matrix
template <typename T>
SparseMatrix<T> csrFileReader<T>::getMatrix(int index) const
{
    if (index < 0 || index >= sparse_matrices.size())
    {
        ModuleBase::WARNING_QUIT("csrFileReader::getMatrix", "Index out of range");
    }
    return sparse_matrices[index];
}

// function to get matrix using R coordinate
template <typename T>
SparseMatrix<T> csrFileReader<T>::getMatrix(int Rx, int Ry, int Rz)
{
    for (int i = 0; i < RCoordinates.size(); i++)
    {
        if (RCoordinates[i][0] == Rx && RCoordinates[i][1] == Ry && RCoordinates[i][2] == Rz)
        {
            return sparse_matrices[i];
        }
    }
    ModuleBase::WARNING_QUIT("csrFileReader::getMatrix", "R coordinate not found");
}

// function to get matrix
template <typename T>
int csrFileReader<T>::getNumberOfR() const
{
    return numberOfR;
}

// function to get matrixDimension
template <typename T>
int csrFileReader<T>::getMatrixDimension() const
{
    return matrixDimension;
}

// function to get step
template <typename T>
int csrFileReader<T>::getStep() const
{
    return step;
}

// T of AtomPair can be double
template class csrFileReader<double>;
// ToDo: T of AtomPair can be std::complex<double>
// template class csrFileReader<std::complex<double>>;

} // namespace ModuleIO