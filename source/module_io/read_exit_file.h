#ifndef READ_EXIT_FILE_H
#define READ_EXIT_FILE_H

#include <string>

namespace ModuleIO
{
/**
 * @brief read file to determine whether to stop the current execution
 *
 * @param my_rank the rank of the current process
 * @param filename the name of the file to read
 * @param ofs_running the output stream
 *
 * @return 0 if the file is not found or does not contain correct keywords, do not stop the current execution
 * @return 1 if the file is found and contains "stop_ion", the current execution stops at the next ionic step
 * @return 2 if the file is found and contains "stop_elec", the current execution stops at the next electronic step
 */
int read_exit_file(const int& my_rank, const std::string& filename, std::ofstream& ofs_running);
} // namespace ModuleIO

#endif