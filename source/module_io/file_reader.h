#ifndef FILE_READER_H
#define FILE_READER_H

#include <fstream>
#include <sstream>
#include <string>

namespace ModuleIO
{

/**
 * @brief A base class of file reader
 * @details This class is supposed to be a base class to read a text file.
 *  it will open a file with a given filename. The function readLine() 
 *  will read a line to a string stream. The function isOpen() check if the file
 *  is open. The destructor will close the file automatically.
 */
class FileReader
{
  public:
    // Default constructor
    FileReader(std::string filename);
    ~FileReader();

    // Check if file is open
    bool isOpen() const;

    // read a line to string stream
    void readLine();

    std::stringstream ss;

  private:
    std::ifstream ifs;
};

} // namespace ModuleIO

#endif