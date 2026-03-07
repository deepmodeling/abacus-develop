#include "restart.h"

#include <fcntl.h>
#include <unistd.h>

#include <fstream>

#include "source_base/global_function.h"
#include "source_base/tool_quit.h"

void Restart::write_file1(const std::string &file_name, const void*const ptr, const size_t size) const
{
	std::ofstream ofs(file_name, std::ofstream::binary|std::ofstream::trunc);
	ofs.write(static_cast<const char*>(ptr),size);
}

void Restart::read_file1(const std::string &file_name, void*const ptr, const size_t size) const
{
	std::ifstream ifs(file_name, std::ifstream::binary);
	ifs.read(static_cast<char*>(ptr),size);
}

bool Restart::write_file2(const std::string& file_name, const void* const ptr, const size_t size, const bool error_quit) const
{
	const int file = open(file_name.c_str(), O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
    if (-1 == file){
        if (error_quit){
            ModuleBase::WARNING_QUIT("Restart::write_file2", "can't open restart save file. errno=" + ModuleBase::GlobalFunc::TO_STRING(errno));
        } else {
            return false;
        }
    }
    auto error = write(file, ptr, size);
    if (-1 == error) {
        if (error_quit) {
            ModuleBase::WARNING_QUIT("Restart::write_file2", "can't write restart save file. errno=" + ModuleBase::GlobalFunc::TO_STRING(errno));
        } else {
            return false;
        }
    }
    close(file);
    return true;
}

namespace GlobalC
{
Restart restart; // Peize Lin add 2020.04.04
} // namespace GlobalC

bool Restart::read_file2(const std::string& file_name, void* const ptr, const size_t size, const bool error_quit) const
{
	const int file = open(file_name.c_str(), O_RDONLY);
    if (-1 == file) {
        if (error_quit) {
            ModuleBase::WARNING_QUIT("Restart::read_file2", "can't open restart load file. errno=" + ModuleBase::GlobalFunc::TO_STRING(errno));
        } else {
            return false;
        }
    }
    auto error = read(file, ptr, size);
    if (-1 == error) {
        if (error_quit) {
            ModuleBase::WARNING_QUIT("Restart::read_file2", "can't read restart load file. errno=" + ModuleBase::GlobalFunc::TO_STRING(errno));
        } else {
            return false;
        }
    }
    close(file);
    return true;
}
