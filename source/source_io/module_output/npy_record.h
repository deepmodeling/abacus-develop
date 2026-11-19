#ifndef ABACUS_NPY_RECORD_H
#define ABACUS_NPY_RECORD_H

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>

namespace ModuleIO
{
// Opens a DOS array and writes its NPY v1 header. Each row starts with the
// returned payload length, schema version and one-based ionic step.
std::array<char, 16> open_npy_record(std::ofstream& out, const std::string& path, bool append, int step,
                                     const std::string& fields, const char* field, std::uint64_t count,
                                     std::uint64_t row_bytes);
}
#endif
