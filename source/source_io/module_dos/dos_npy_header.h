#ifndef MODULE_IO_DOS_NPY_HEADER_H
#define MODULE_IO_DOS_NPY_HEADER_H

#include <algorithm>
#include <cstdint>
#include <limits>
#include <ostream>
#include <string>

namespace ModuleIO
{
/// Write the NPY v1 header for a one-dimensional array of structured records.
inline void write_dos_npy_header(std::ostream& out, const std::string& fields, std::size_t count)
{
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "DOS requires binary64");
    const std::uint16_t one = 1;
    const char endian = *reinterpret_cast<const unsigned char*>(&one) ? '<' : '>';
    std::string header = "{'descr': [" + fields + "], 'fortran_order': False, 'shape': ("
                       + std::to_string(count) + ",), }";
    std::replace(header.begin(), header.end(), '=', endian);
    header.append((64 - (11 + header.size()) % 64) % 64, ' ');
    header += '\n';
    if (header.size() > std::numeric_limits<std::uint16_t>::max())
    {
        out.setstate(std::ios::failbit);
        return;
    }
    const unsigned char length[] = {static_cast<unsigned char>(header.size()),
                                    static_cast<unsigned char>(header.size() >> 8)};
    out.write("\x93NUMPY\x01\x00", 8);
    out.write(reinterpret_cast<const char*>(length), sizeof(length));
    out.write(header.data(), header.size());
}
}

#endif
