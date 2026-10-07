#ifndef POS_OP_CSR_H
#define POS_OP_CSR_H

#include <string>

namespace ModuleIO
{
namespace detail
{
bool lat_r_nonempty(const int nonzero_num[3]);

void assemble_csr(const std::string& output_filename,
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
