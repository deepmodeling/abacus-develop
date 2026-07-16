#include "rr_sparse_writer.h"

#include "source_base/tool_quit.h"

#include <fstream>

namespace ModuleIO
{
namespace detail
{
bool rr_sparse_has_payload(const int nonzero_num[3])
{
    return nonzero_num[0] != 0 || nonzero_num[1] != 0 || nonzero_num[2] != 0;
}

void finalize_rr_sparse_file(const std::string& output_filename,
                             const std::string& payload_filename,
                             const int step,
                             const int nlocal,
                             const int output_R_number,
                             const bool binary,
                             const bool append,
                             const std::string& context)
{
    std::ios_base::openmode output_mode = std::ios::out;
    std::ios_base::openmode payload_mode = std::ios::in;
    if (binary)
    {
        output_mode |= std::ios::binary;
        payload_mode |= std::ios::binary;
    }
    if (append)
    {
        output_mode |= std::ios::app;
    }

    std::ofstream out_r(output_filename.c_str(), output_mode);
    if (!out_r.is_open())
    {
        ModuleBase::WARNING_QUIT(context, "Cannot open r(R) output file: " + output_filename);
    }

    if (binary)
    {
        out_r.write(reinterpret_cast<const char*>(&step), sizeof(int));
        out_r.write(reinterpret_cast<const char*>(&nlocal), sizeof(int));
        out_r.write(reinterpret_cast<const char*>(&output_R_number), sizeof(int));
    }
    else
    {
        out_r << "STEP: " << step << std::endl;
        out_r << "Matrix Dimension of r(R): " << nlocal << std::endl;
        out_r << "Matrix number of r(R): " << output_R_number << std::endl;
    }

    std::ifstream payload(payload_filename.c_str(), payload_mode);
    if (!payload.is_open())
    {
        ModuleBase::WARNING_QUIT(context, "Cannot read temporary sparse matrix file: " + payload_filename);
    }
    out_r << payload.rdbuf();
}
} // namespace detail
} // namespace ModuleIO
