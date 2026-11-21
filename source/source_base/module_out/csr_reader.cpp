#include "csr_reader.h"
#include <cmath>
#include <set>
#include <limits>
#include "source_base/tool_quit.h"

namespace ModuleIO
{

namespace
{
/**
 * @brief Extract the ionic step from a frame header line held in `ss`.
 * @details Accepts the current bare "Ionic_Step N" spelling, the legacy
 * compact "STEP: N" label and the historical "--- Ionic Step N ---" /
 * "--- Ionic_Step N ---" forms, so files written by older versions remain
 * readable.
 * @param ss stream positioned at the start of the header line
 * @param step output ionic step index
 * @return true when a recognized marker and step value were parsed
 */
bool parseStepMarker(std::stringstream& ss, int& step)
{
    std::string first;
    if (!(ss >> first))
    {
        return false;
    }
    if (first == "STEP:" || first == "Ionic_Step")
    {
        return static_cast<bool>(ss >> step);
    }
    if (first != "---")
    {
        return false;
    }
    std::string ionic;
    if (!(ss >> ionic))
    {
        return false;
    }
    if (ionic == "Ionic_Step")
    {
        return static_cast<bool>(ss >> step);
    }
    if (ionic == "Ionic")
    {
        std::string label;
        return static_cast<bool>(ss >> label >> step) && label == "Step";
    }
    return false;
}
} // namespace

// constructor
template <typename T>
csrFileReader<T>::csrFileReader(const std::string& filename) : FileReader(filename)
{
    parseFile();
}

// Format-specific metadata is followed by the same CSR block reader.
template <typename T>
void csrFileReader<T>::parseFile()
{
    RCoordinates.clear();
    sparse_matrices.clear();
    nextDataLine();
    const std::string line = ss.str();
    const std::size_t begin = line.find_first_not_of(" \t");
    bool compact = begin != std::string::npos && line.compare(begin, 5, "STEP:") == 0;
    if (!compact)
    {
        // The bare "Ionic_Step N" marker is shared by the compact and
        // cell-aware headers; peek at the next physical line to tell them
        // apart. Compact headers continue with "Matrix Dimension of ...",
        // while cell-aware headers continue with a description comment.
        std::string next;
        if (peekLine(next))
        {
            const std::size_t nbegin = next.find_first_not_of(" \t");
            compact = nbegin != std::string::npos
                      && next.compare(nbegin, 19, "Matrix Dimension of") == 0;
        }
    }
    if (compact)
    {
        parseCompactHeader();
    }
    else
    {
        parseCellHeader();
    }
    if (matrixDimension <= 0 || numberOfR < 0)
    {
        ModuleBase::WARNING_QUIT("csrFileReader::parseFile", "Invalid CSR dimensions or R-block count");
    }
    readBlocks();
    validateFrameEnd();
}

template <typename T>
void csrFileReader<T>::parseCompactHeader()
{
    const auto invalid = []() {
        ModuleBase::WARNING_QUIT("csrFileReader::parseCompactHeader", "Malformed compact CSR file");
    };
    if (!parseStepMarker(ss, step)) { invalid(); }
    const auto header_integer = [&](const std::string& prefix, int& value) {
        readLine();
        const std::string line = ss.str();
        const std::size_t begin = line.find_first_not_of(" \t");
        const std::size_t colon = line.find(':');
        if (begin == std::string::npos || line.compare(begin, prefix.size(), prefix) != 0
            || colon == std::string::npos) { invalid(); }
        std::istringstream values(line.substr(colon + 1));
        if (!(values >> value)) { invalid(); }
    };
    header_integer("Matrix Dimension of ", matrixDimension);
    header_integer("Matrix number of ", numberOfR);
}

// Cell-aware headers used by HContainer writers include spin and geometry.
template <typename T>
void csrFileReader<T>::parseCellHeader()
{
    const auto invalid = []() {
        ModuleBase::WARNING_QUIT("csrFileReader::parseCellHeader", "Malformed cell-aware CSR header");
    };
    if (!parseStepMarker(ss, step)) { invalid(); }
    readLine(); // matrix description
    readLine();
    int spins = 0;
    if (!(ss >> spins) || spins <= 0) { invalid(); }
    readLine();
    int spin = 0;
    if (!(ss >> spin) || spin < 1 || spin > spins) { invalid(); }
    readLine();
    if (!(ss >> matrixDimension)) { invalid(); }
    readLine();
    if (!(ss >> numberOfR)) { invalid(); }
    nextDataLine(); // lattice name, possibly after an empty line
    for (int i = 0; i < 5; ++i) { readLine(); } // scale, three vectors, species
    readLine(); // atom counts
    int atoms = 0;
    int count = 0;
    while (true)
    {
        ss >> std::ws;
        if (ss.eof() || ss.peek() == '#') { break; }
        if (!(ss >> count) || count < 0 || count > std::numeric_limits<int>::max() - atoms)
        { invalid(); }
        atoms += count;
    }
    if (atoms <= 0) { invalid(); }
    readLine(); // coordinate convention
    for (int i = 0; i < atoms; ++i) { readLine(); }
}

// Skip only whitespace and comments, never malformed numerical data.
template <typename T>
void csrFileReader<T>::nextDataLine()
{
    do
    {
        readLine();
        ss >> std::ws;
    } while (ss.eof() || ss.peek() == '#');
}

template <typename T>
void csrFileReader<T>::nextDataToken()
{
    ss >> std::ws;
    if (ss.eof() || ss.peek() == '#') { nextDataLine(); }
}

template <typename T>
void csrFileReader<T>::readBlocks()
{
    const auto invalid = []() {
        ModuleBase::WARNING_QUIT("csrFileReader::readBlocks", "Malformed CSR block");
    };
    std::set<std::vector<int>> seen;
    for (int ir = 0; ir < numberOfR; ++ir)
    {
        nextDataLine();
        std::vector<int> coordinate(3);
        int nnz = 0;
        if (!(ss >> coordinate[0] >> coordinate[1] >> coordinate[2] >> nnz) || nnz < 0)
        { invalid(); }
        const std::size_t dimension = matrixDimension;
        if (static_cast<std::size_t>(nnz) > dimension * dimension || !seen.insert(coordinate).second)
        { invalid(); }
        std::vector<T> values(nnz);
        std::vector<int> columns(nnz);
        std::vector<int> pointers(dimension + 1);
        for (int i = 0; i < nnz; ++i)
        {
            nextDataToken();
            if (!(ss >> values[i]) || !std::isfinite(std::abs(values[i]))) { invalid(); }
        }
        for (int i = 0; i < nnz; ++i)
        {
            nextDataToken();
            if (!(ss >> columns[i]) || columns[i] < 0 || columns[i] >= matrixDimension) { invalid(); }
        }
        for (int i = 0; i <= matrixDimension; ++i)
        {
            nextDataToken();
            if (!(ss >> pointers[i]) || pointers[i] < 0 || pointers[i] > nnz) { invalid(); }
            if (i > 0 && pointers[i] < pointers[i - 1]) { invalid(); }
        }
        if (pointers.front() != 0 || pointers.back() != nnz) { invalid(); }
        ss >> std::ws;
        if (!ss.eof() && ss.peek() != '#') { invalid(); }
        for (int row = 0; row < matrixDimension; ++row)
        {
            for (int i = pointers[row] + 1; i < pointers[row + 1]; ++i)
            {
                if (columns[i] <= columns[i - 1]) { invalid(); }
            }
        }
        SparseMatrix<T> matrix(matrixDimension, matrixDimension);
        matrix.readCSR(values, columns, pointers);
        RCoordinates.push_back(coordinate);
        sparse_matrices.push_back(matrix);
    }
}

// Only blank/comment lines, EOF or an unconsumed next frame may follow a frame.
template <typename T>
void csrFileReader<T>::validateFrameEnd()
{
    std::string line;
    while (peekLine(line))
    {
        std::istringstream tokens(line);
        tokens >> std::ws;
        if (tokens.eof() || tokens.peek() == '#')
        {
            readLine();
            continue;
        }
        std::string marker;
        tokens >> marker;
        if (marker == "STEP:" || marker == "Ionic_Step") { return; }
        std::string ionic;
        if (marker == "---" && tokens >> ionic)
        {
            if (ionic == "Ionic_Step") { return; }
            std::string step_label;
            if (ionic == "Ionic" && tokens >> step_label && step_label == "Step") { return; }
        }
        ModuleBase::WARNING_QUIT("csrFileReader::validateFrameEnd", "Unexpected data after CSR frame");
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
