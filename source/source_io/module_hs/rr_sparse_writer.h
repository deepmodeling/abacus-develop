#ifndef RR_SPARSE_WRITER_H
#define RR_SPARSE_WRITER_H

#include <string>

namespace ModuleIO
{
namespace detail
{
bool rr_sparse_has_payload(const int nonzero_num[3]);

void finalize_rr_sparse_file(const std::string& output_filename,
                             const std::string& payload_filename,
                             int step,
                             int nlocal,
                             int output_R_number,
                             bool binary,
                             bool append,
                             const std::string& context);
} // namespace detail
} // namespace ModuleIO

#endif
