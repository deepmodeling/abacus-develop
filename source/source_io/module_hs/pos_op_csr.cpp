#include "pos_op_csr.h"

#include "source_base/tool_quit.h"
#include "source_cell/ucell_io.h"
#include "source_io/module_hs/lat_r_csr.h"
#include "source_io/module_output/provenance_header.h"

#include <fstream>

namespace ModuleIO
{
namespace detail
{
bool lat_r_nonempty(const int nonzero_num[3])
{
    return nonzero_num[0] != 0 || nonzero_num[1] != 0 || nonzero_num[2] != 0;
}

void assemble_csr(const std::string& output_filename,
                  const std::string& payload_filename,
                  const int step,
                  const int nlocal,
                  const int output_R_number,
                  const bool binary,
                  const bool append,
                  const std::string& axis,
                  const UnitCell* ucell,
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
        if (ucell == nullptr)
        {
            ModuleBase::WARNING_QUIT(context, "UnitCell pointer must not be null for text r(R) output.");
        }
        const std::string time_stamp = ModuleIO::local_timestamp();
        out_r << "Ionic_Step " << step << ModuleIO::prov_inline(time_stamp) << std::endl;
        out_r << " # print r matrix in real space r(R) " << axis << "-component" << std::endl;
        out_r << " 1 # number of spin directions" << std::endl;
        out_r << " 1 # spin index" << std::endl;
        out_r << " " << nlocal << " # number of localized basis" << std::endl;
        out_r << " " << output_R_number << " # number of Bravais lattice vector R" << std::endl;
        out_r << std::endl;
        ModuleIO::UcellIO::write_ucell(out_r, ucell);
        out_r << std::endl;
        ModuleIO::write_csr_banner(out_r);
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
