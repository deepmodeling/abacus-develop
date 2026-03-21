#include "csr_reader.h"
#include "source_base/tool_quit.h"

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
    ss >> tmp_string >> tmp_string >> tmp_string >> step;

    //std::cout << " step is " << step << std::endl; 
    // Read the title
    readLine();
    // Read the total spin
    readLine();
    // Read the spin index
    readLine();

    // Read the matrix dimension
    readLine();
    ss >> matrixDimension;
    // std::cout << " mat dim is " << matrixDimension << std::endl; 

    // Read the number of R
    readLine();
    ss >> numberOfR;
    // std::cout << " number of R is " << numberOfR << std::endl;
    readLine();

    // Read cell
    read_ucell();

    // Read CSR format
    readLine();
    readLine();
    readLine();
    readLine();
    readLine();
    readLine();
    readLine();
    readLine();
    readLine(); // read the last line of CSR format

    // Read the matrices
    for (int i = 0; i < numberOfR; i++)
    {
        // std::cout << " read R " << i+1 << std::endl;

        std::vector<int> RCoord(3);
        int nonZero = 0;

        readLine();
        readLine();
        ss >> RCoord[0] >> RCoord[1] >> RCoord[2] >> nonZero;
        RCoordinates.push_back(RCoord);

        std::vector<T> csr_values(nonZero);
        std::vector<int> csr_col_ind(nonZero);
        std::vector<int> csr_row_ptr(matrixDimension + 1);

        // read CSR values
        readLine();
        // std::cout << " ss1: " << ss.str() << std::endl;

        readLine();
	size_t count1 = 0;
        while (count1 < nonZero)
        {
            if (ss.eof() || ss.fail())
            {
                readLine();
	    }
            if (ss >> csr_values[count1])
            {
                count1++;
            }
	}
        // std::cout << "count1=" << count1 << std::endl;

        // read CSR column indices
        readLine();
        // std::cout << " ss2: " << ss.str() << std::endl;

	size_t count2 = 0;
        while (count2 < nonZero)
        {
            if (ss.eof() || ss.fail())
            {
                readLine();
	    }
            if (ss >> csr_col_ind[count2])
            {
                count2++;
            }
	}
        // std::cout << "count2=" << count2 << std::endl;

        // read row pointers
        readLine();
        // std::cout << " ss3: " << ss.str() << std::endl;

	size_t count3 = 0;
        while (count3 < matrixDimension + 1)
        {
            if (ss.eof() || ss.fail())
            {
                readLine();
	    }
            if (ss >> csr_row_ptr[count3])
            {
                count3++;
            }
	}
        // std::cout << "count3=" << count3 << std::endl;

        // create sparse matrix
        SparseMatrix<T> matrix(matrixDimension, matrixDimension);
        matrix.readCSR(csr_values, csr_col_ind, csr_row_ptr);
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
SparseMatrix<T> csrFileReader<T>::getMatrix(int Rx, int Ry, int Rz) const
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
template class csrFileReader<std::complex<double>>;

} // namespace ModuleIO
