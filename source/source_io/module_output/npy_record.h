#ifndef ABACUS_NPY_RECORD_H
#define ABACUS_NPY_RECORD_H

#include <array>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <string>

namespace ModuleIO
{
// Every row starts with the frame's payload size (little endian), schema
// version and one-based ionic step. Remaining numeric fields use native order.
using NpyPrefix = std::array<char, 16>;
using NpyRows = std::function<void(std::ostream&, const NpyPrefix&)>;

// fields is a comma-separated Python dtype field list. '=' in its numeric
// descriptors is replaced by the writer's explicit byte order.
void write_npy_record(const std::string& path, bool append, int step,
                      const std::string& fields, const std::string& identifying_field,
                      std::uint64_t count, std::uint64_t row_bytes, const NpyRows& write_rows);
}
#endif
