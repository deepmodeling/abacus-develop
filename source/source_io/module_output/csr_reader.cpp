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

    // Read cell
    // Skip empty line before ucell info if present
    readLine();
    std::string line = ss.str();
    if (line.empty() || line.find_first_not_of(" \t\n\r") == std::string::npos)
    {
        // Empty line, read next line for latName
        readLine();
    }
    // Now ss contains latName (or the line after empty line)
    // We don't need to use latName, just continue reading
    
    // Read lat0, latvec (3 lines), atom labels (6 lines total including latName)
    readLine(); // lat0
    readLine(); // latvec e1
    readLine(); // latvec e2
    readLine(); // latvec e3
    readLine(); // atom labels
    
    // Read atom numbers
    readLine();
    std::istringstream iss(ss.str());
    int natom = 0;
    int total_natom = 0;
    while (iss >> natom)
    {
        total_natom += natom;
    }
    
    // Read "Direct" line and atom coordinates
    for (int i = 0; i < total_natom + 1; i++)
    {
        readLine(); // Direct + atom coordinates
    }

    // Skip empty lines and CSR format comment block
    // (lines starting with '#' or containing only whitespace)
    // Use readLine() since ifs is private in FileReader
    bool found_data = false;
    while (!found_data)
    {
        readLine();
        std::string line = ss.str();
        size_t start = line.find_first_not_of(" \t");
        if (start != std::string::npos && line[start] != '#')
        {
            // This is a data line (first R-block), ss already holds it
            found_data = true;
        }
    }

    // Read the matrices
    // Note: ss already contains the first R-block line from the skip loop above
    for (int i = 0; i < numberOfR; i++)
    {
        // std::cout << " read R " << i+1 << std::endl;

        std::vector<int> RCoord(3);
        int nonZero = 0;

        if (i > 0)
        {
            // For subsequent R-blocks, skip empty lines to find next R-coordinate line
            bool found = false;
            while (!found)
            {
                readLine();
                std::string line = ss.str();
                size_t start = line.find_first_not_of(" \t");
                if (start != std::string::npos && line[start] != '#')
                {
                    found = true;
                }
            }
        }
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
