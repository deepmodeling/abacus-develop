#include "csr_reader.h"
#include "source_base/tool_quit.h"
#include <cctype>
#include <sstream>
#include <string>

namespace
{
bool is_digit_char(const char c)
{
    return std::isdigit(static_cast<unsigned char>(c)) != 0;
}

int parse_last_int_token(const std::string& line)
{
    std::stringstream iss(line);
    std::string token;
    int value = 0;
    bool found = false;
    while (iss >> token)
    {
        size_t begin = 0;
        while (begin < token.size() && !(is_digit_char(token[begin]) || token[begin] == '-' || token[begin] == '+'))
        {
            ++begin;
        }
        if (begin == token.size())
        {
            continue;
        }
        size_t end = token.size();
        while (end > begin && !is_digit_char(token[end - 1]))
        {
            --end;
        }
        try
        {
            value = std::stoi(token.substr(begin, end - begin));
            found = true;
        }
        catch (...)
        {
        }
    }
    if (!found)
    {
        ModuleBase::WARNING_QUIT("csrFileReader::parseFile", "Failed to parse integer from line: " + line);
    }
    return value;
}
}

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

    std::string tmp_string;

    readLine();
    step = parse_last_integer(ss.str(), "csrFileReader::parseFile");

    // Read the step and detect file flavor.
    readLine();
    const std::string first_line = ss.str();
    const bool compact_format = first_line.rfind("STEP:", 0) == 0;
    if (compact_format)
    {
        step = parse_last_int_token(first_line);

        readLine();
        matrixDimension = parse_last_int_token(ss.str());

        readLine();
        numberOfR = parse_last_int_token(ss.str());
    }
    else
    {
        ss >> tmp_string >> tmp_string >> tmp_string >> step;

        // Read the title
        readLine();
        // Read the total spin
        readLine();
        // Read the spin index
        readLine();

        // Read the matrix dimension
        readLine();
        ss >> matrixDimension;

        // Read the number of R
        readLine();
        ss >> numberOfR;
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
    }

    // Read the matrices
    // Note: ss already contains the first R-block line from the skip loop above
    for (int i = 0; i < numberOfR; i++)
    {
        std::vector<int> RCoord(3);
        int nonZero = 0;

        if (!compact_format)
        {
            readLine();
        }
        readLine();
        ss >> RCoord[0] >> RCoord[1] >> RCoord[2] >> nonZero;
        RCoordinates.push_back(RCoord);

        std::vector<T> csr_values(nonZero);
        std::vector<int> csr_col_ind(nonZero);
        std::vector<int> csr_row_ptr(matrixDimension + 1);

        // read CSR values
        if (!compact_format)
        {
            readLine();
        }
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
        if (!compact_format)
        {
            readLine();
        }
        readLine();

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
        if (!compact_format)
        {
            readLine();
        }
        readLine();

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
