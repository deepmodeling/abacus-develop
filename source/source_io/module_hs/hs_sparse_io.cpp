#include "hs_sparse_io.h"
#include "hs_sparse_io_detail.h"
#include "source_io/module_output/provenance_header.h"

#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/ucell_io.h"
#include "lat_r_csr.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

namespace
{
int count_output_R(const std::vector<long long>& nonzero_num)
{
    int output_R_number = 0;
    for (const long long count: nonzero_num)
    {
        if (count != 0)
        {
            ++output_R_number;
        }
    }
    return output_R_number;
}

void open_sparse_file(std::ofstream& ofs, const ModuleIO::SparseWriteOptions& options)
{
    std::ios_base::openmode mode = std::ios::out;
    if (options.binary)
    {
        mode |= std::ios::binary;
    }
    // Append from the second ionic step onward in append mode; the first
    // step (istep == 0) truncates so a rerun does not retain prior frames.
    // Single-point SCF passes istep = -1 and also starts a fresh file.
    const bool append_on_restart = (options.out_app_flag == 1) && (options.istep > 0);
    if (append_on_restart)
    {
        mode |= std::ios::app;
    }
    ofs.open(options.filename.c_str(), mode);
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::open_sparse_file",
                                 "Cannot open sparse matrix file: " + options.filename);
    }
}

void write_sparse_header(std::ofstream& ofs,
                         const ModuleIO::SparseWriteOptions& options,
                         const int nlocal,
                         const int output_R_number)
{
    // ionic step printed one-based per governance rule 22, consistent
    // with the g{N} filename suffix (g1 holds the first ionic step)
    const int step = std::max(options.istep, 0) + 1;
    if (options.binary)
    {
        ofs.write(reinterpret_cast<const char*>(&step), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&nlocal), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&output_R_number), sizeof(int));
    }
    else
    {
        const std::string time_stamp = ModuleIO::local_timestamp();
        ofs << "Ionic_Step " << step << ModuleIO::prov_inline(time_stamp) << std::endl;
        ofs << "Matrix Dimension of " + options.label + "(R): " << nlocal
            << std::endl;
        ofs << "Matrix number of " + options.label + "(R): "
            << output_R_number << std::endl;
    }
}

void write_R_record(std::ofstream& ofs,
                    const ModuleIO::RCoordinate& R_coor,
                    const long long nonzero_count,
                    const bool binary,
                    const bool leading_space)
{
    int dRx = R_coor.x;
    int dRy = R_coor.y;
    int dRz = R_coor.z;
    if (binary)
    {
        const int count = static_cast<int>(nonzero_count);
        ofs.write(reinterpret_cast<char*>(&dRx), sizeof(int));
        ofs.write(reinterpret_cast<char*>(&dRy), sizeof(int));
        ofs.write(reinterpret_cast<char*>(&dRz), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&count), sizeof(int));
    }
    else
    {
        if (leading_space)
        {
            ofs << " ";
        }
        ofs << dRx << " " << dRy << " " << dRz << " " << nonzero_count
            << std::endl;
    }
}

/**
 * @brief Write the four-line record used by the H(R)/S(R) writer for an R
 *        block whose elements are all below the sparse threshold.
 */
void write_empty_csr_block(std::ofstream& ofs, const int nbasis)
{
    ofs << " # CSR values" << std::endl;
    ofs << std::endl;
    ofs << " # CSR column indices" << std::endl;
    ofs << std::endl;
    ofs << " # CSR row pointers" << std::endl;
    for (int i = 0; i < nbasis + 1; ++i)
    {
        ofs << " 0";
    }
    ofs << std::endl;
}

void write_csr_style_header(std::ofstream& ofs,
                            const ModuleIO::SparseWriteOptions& options,
                            const int nlocal,
                            const int output_R_number,
                            const UnitCell* ucell)
{
    if (ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_sparse_csr",
                                 "UnitCell pointer must not be null for text output.");
    }
    // ionic step printed one-based per governance rule 22, consistent
    // with the g{N} filename suffix (g1 holds the first ionic step)
    const int step = std::max(options.istep, 0) + 1;
    const std::string time_stamp = ModuleIO::local_timestamp();
    ofs << "Ionic_Step " << step << ModuleIO::prov_inline(time_stamp) << std::endl;
    ofs << " # print " << options.label << " matrix in real space "
        << options.label << "(R)" << std::endl;
    ofs << " " << options.nspin << " # number of spin directions" << std::endl;
    ofs << " " << options.ispin + 1 << " # spin index" << std::endl;
    ofs << " " << nlocal << " # number of localized basis" << std::endl;
    ofs << " " << output_R_number << " # number of Bravais lattice vector R" << std::endl;
    ofs << std::endl;
    ModuleIO::UcellIO::write_ucell(ofs, ucell);
    ofs << std::endl;
    ModuleIO::write_csr_banner(ofs);
}
} // namespace

template <typename Tdata>
void ModuleIO::save_sparse(
    const SparseRMatrix<Tdata>& smat,
    const std::set<RCoordinate>& all_R_coor,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options) {
    ModuleBase::TITLE("ModuleIO", "save_sparse");
    ModuleBase::timer::start("ModuleIO", "save_sparse");
    const int nlocal = pv.para_2d().nrglb();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_sparse",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    const std::vector<long long> nonzero_num
        = detail::count_nonzeros_by_R(smat, all_R_coor, options.threshold, options.reduce);
    const bool text_csr = options.csr_format && !options.binary;
    // The H(R)/S(R) CSR files list every R vector and write an all-zero block
    // for empty ones, so the header count and the per-block records must cover
    // the full R set. The compact legacy layout still skips empty blocks.
    const int output_R_number = text_csr ? static_cast<int>(all_R_coor.size())
                                         : count_output_R(nonzero_num);
    std::ofstream ofs;
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        open_sparse_file(ofs, options);
        if (text_csr)
        {
            write_csr_style_header(ofs, options, nlocal, output_R_number, options.ucell);
        }
        else
        {
            write_sparse_header(ofs, options, nlocal, output_R_number);
        }
    }

    const SparseRBlock<Tdata> empty_block;
    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        const bool is_empty = (nonzero_num[count] == 0);
        // The compact legacy layout skips empty R blocks entirely. The CSR
        // text layout lists every R vector and writes a zero block for empty
        // ones, matching the H(R)/S(R) writer.
        if (!text_csr && is_empty)
        {
            count++;
            continue;
        }

        if (!options.reduce || GlobalV::DRANK == 0)
        {
            write_R_record(ofs, R_coor, nonzero_num[count], options.binary, text_csr);
            if (text_csr && is_empty)
            {
                write_empty_csr_block(ofs, nlocal);
            }
        }

        if (text_csr && is_empty)
        {
            ++count;
            continue;
        }

        const SparseRBlock<Tdata>& block = (smat.count(R_coor) != 0) ? smat.at(R_coor) : empty_block;
        if (text_csr)
        {
            save_lat_r_csr(ofs, block, pv.para_2d(), options);
        }
        else
        {
            save_lat_r(ofs, block, pv.para_2d(), options);
        }
        ++count;
    }
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        ofs.close();
    }

    ModuleBase::timer::end("ModuleIO", "save_sparse");
}

template void ModuleIO::save_sparse<double>(
    const SparseRMatrix<double>&,
    const std::set<RCoordinate>&,
    const Parallel_Orbitals&,
    const SparseWriteOptions&);

template void ModuleIO::save_sparse<std::complex<double>>(
    const SparseRMatrix<std::complex<double>>&,
    const std::set<RCoordinate>&,
    const Parallel_Orbitals&,
    const SparseWriteOptions&);
