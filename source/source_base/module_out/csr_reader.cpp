#include "csr_reader.h"
#include "source_base/tool_quit.h"

namespace ModuleIO
{
namespace
{
bool is_blank_line(const std::string& line)
{
    return line.find_first_not_of(" \t\r") == std::string::npos;
}

bool is_comment_line(const std::string& line)
{
    const std::size_t first = line.find_first_not_of(" \t\r");
    return first != std::string::npos && line[first] == '#';
}

bool contains_token(const std::string& line, const std::string& token)
{
    return line.find(token) != std::string::npos;
}

int parse_last_integer(const std::string& line, const std::string& context)
{
    std::stringstream parser(line);
    int value = 0;
    bool found = false;
    while (!parser.eof())
    {
        if (parser >> value)
        {
            found = true;
            continue;
        }
        parser.clear();
        std::string dummy;
        parser >> dummy;
    }
    if (!found)
    {
        ModuleBase::WARNING_QUIT(context, "Failed to parse integer from line: " + line);
    }
    return value;
}

void read_next_payload_line(FileReader& reader)
{
    do
    {
        reader.readLine();
    } while (is_blank_line(reader.ss.str()) || is_comment_line(reader.ss.str()));
}

template <typename Tvalue>
void read_numeric_block(FileReader& reader, std::vector<Tvalue>& buffer)
{
    std::size_t count = 0;
    while (count < buffer.size())
    {
        read_next_payload_line(reader);
        while (count < buffer.size() && (reader.ss >> buffer[count]))
        {
            ++count;
        }
    }
}
} // namespace


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

    readLine();
    step = parse_last_integer(ss.str(), "csrFileReader::parseFile");

    // Support both the legacy verbose CSR header and the compact
    // Matrix-Dimension/Matrix-number header produced by save_sparse().
    readLine();
    const std::string second_line = ss.str();
    if (contains_token(second_line, "Matrix Dimension of"))
    {
        matrixDimension = parse_last_integer(second_line, "csrFileReader::parseFile");
        readLine();
        numberOfR = parse_last_integer(ss.str(), "csrFileReader::parseFile");
    }
    else
    {
        // Legacy verbose header:
        // title -> total spin -> spin index -> matrix dimension -> number of R
        readLine();
        readLine();
        readLine();
        matrixDimension = parse_last_integer(ss.str(), "csrFileReader::parseFile");
        readLine();
        numberOfR = parse_last_integer(ss.str(), "csrFileReader::parseFile");
        readLine();
        read_ucell();
    }

    // Read the matrices
    // Note: ss already contains the first R-block line from the skip loop above
    for (int i = 0; i < numberOfR; i++)
    {
        std::vector<int> RCoord(3);
        int nonZero = 0;

        read_next_payload_line(*this);
        ss >> RCoord[0] >> RCoord[1] >> RCoord[2] >> nonZero;
        RCoordinates.push_back(RCoord);

        std::vector<T> csr_values(nonZero);
        std::vector<int> csr_col_ind(nonZero);
        std::vector<int> csr_row_ptr(matrixDimension + 1);

        read_numeric_block(*this, csr_values);
        read_numeric_block(*this, csr_col_ind);
        read_numeric_block(*this, csr_row_ptr);

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
