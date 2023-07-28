#include "file_reader.h"

#include "module_base/tool_quit.h"

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

} // namespace ModuleIO