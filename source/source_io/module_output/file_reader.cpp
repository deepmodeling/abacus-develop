#include "file_reader.h"

#include "source_base/tool_quit.h"

namespace ModuleIO
{

// Constructor
FileReader::FileReader(std::string filename)
{
    ifs.open(filename.c_str());
    if (!ifs.is_open())
    {
        ModuleBase::WARNING_QUIT("FileReader::FileReader", "Error opening file");
    }
}

// Destructor
FileReader::~FileReader()
{
    if (ifs.is_open())
    {
        ifs.close();
    }
}

// Function to check if file is open
bool FileReader::isOpen() const
{
    return ifs.is_open();
}

// Function to read a line and return string stream
void FileReader::readLine()
{
    // add warning if file is not open
    if (!ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        ss.clear();
        ss.str(line);
    }
    else
    {
        ModuleBase::WARNING_QUIT("FileReader::readLine", "End of file");
    }
}

void FileReader::read_ucell()
{
    if (ifs.eof())
    {
        ModuleBase::WARNING_QUIT("FileReader::read_ucell", "End of file");
    }

    std::string tmp;
    for (int i = 0; i < 6; i++)
    {
        std::getline(ifs, tmp); // latName + lat0 + latvec + atom label
    }
    std::getline(ifs, tmp); // atom number of each type

    std::istringstream iss(tmp);
    int natom = 0;
    int total_natom = 0;
    while (iss >> natom)
    {
        total_natom += natom;
    }
    for (int i = 0; i < total_natom + 1; i++)
    {
        std::getline(ifs, tmp); // Direct + atom coordinates
    }
}


} // namespace ModuleIO
