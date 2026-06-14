#include "restart.h"

#include <fcntl.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <io.h> // open/read/write/close (_open ...) on Windows
// The POSIX owner-permission bits are not defined by the Windows CRT; map them
// to the MSVCRT read/write mode bits so open(..., O_CREAT, mode) still works.
#ifndef S_IRUSR
#define S_IRUSR _S_IREAD
#endif
#ifndef S_IWUSR
#define S_IWUSR _S_IWRITE
#endif
#else
#include <unistd.h>
#endif

#include <fstream>
#include <stdexcept>

#include "source_base/global_function.h"

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
            throw std::runtime_error("can't open restart save file. \nerrno=" + ModuleBase::GlobalFunc::TO_STRING(errno) + ".\n" + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        } else {
            return false;
        }
    }
    auto error = write(file, ptr, size);
    if (-1 == error) {
        if (error_quit) {
            throw std::runtime_error("can't write restart save file. \nerrno=" + ModuleBase::GlobalFunc::TO_STRING(errno) + ".\n" + std::string(__FILE__) + " line " + std::to_string(__LINE__));
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
            throw std::runtime_error("can't open restart load file. \nerrno=" + ModuleBase::GlobalFunc::TO_STRING(errno) + ".\n" + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        } else {
            return false;
        }
    }
    auto error = read(file, ptr, size);
    if (-1 == error) {
        if (error_quit) {
            throw std::runtime_error("can't read restart load file. \nerrno=" + ModuleBase::GlobalFunc::TO_STRING(errno) + ".\n" + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        } else {
            return false;
        }
    }
    close(file);
    return true;
}