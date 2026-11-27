#include "hs_sparse_io.h"
#include "hs_sparse_io_detail.h"
#include "source_io/module_output/provenance_header.h"

#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "lat_r_csr.h"

#include <algorithm>
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
    const bool append_on_restart = ModuleIO::detail::append_ionic_step(options.out_app_flag,
                                                                      options.istep);
    if (append_on_restart)
    {
        mode |= std::ios::app;
    }
    ofs.open(options.filename.c_str(), mode);
    ModuleIO::detail::check_output_file_open(ofs, options.filename, "ModuleIO::open_sparse_file");
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
        ModuleIO::detail::write_csr_binary_header(ofs, step, nlocal, output_R_number);
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

void write_csr_style_header(std::ofstream& ofs,
                            const ModuleIO::SparseWriteOptions& options,
                            const int nlocal,
                            const int output_R_number)
{
    // ionic step printed one-based per governance rule 22, consistent
    // with the g{N} filename suffix (g1 holds the first ionic step)
    const int step = std::max(options.istep, 0) + 1;
    const std::string description
        = " # print " + options.label + " matrix in real space " + options.label + "(R)";
    const int ispin_one_based = options.ispin + 1;

    const ModuleIO::detail::CsrTextHeader header = {
        step, description, options.nspin, ispin_one_based, nlocal, output_R_number,
        options.ucell, std::string(), false, 0.0, true};
    ModuleIO::detail::write_csr_text_header(ofs, header);
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
    // Sparse output only lists R blocks that carry at least one element above
    // the threshold; the header count therefore matches the number of blocks
    // actually written for both the self-describing CSR layout and the compact
    // legacy layout.
    const int output_R_number = count_output_R(nonzero_num);
    std::ofstream ofs;
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        open_sparse_file(ofs, options);
        if (text_csr)
        {
            write_csr_style_header(ofs, options, nlocal, output_R_number);
        }
        else
        {
            write_sparse_header(ofs, options, nlocal, output_R_number);
        }
    }

    const SparseRBlock<Tdata> empty_block;
    int count = 0;
    for (const Abfs::Vector3_Order<int>& R_coor: all_R_coor)
    {
        const bool is_empty = (nonzero_num[count] == 0);
        // Sparse output must not persist all-zero blocks: they carry no
        // information and make files grow with the candidate R grid rather
        // than with the non-zero R set.
        if (is_empty)
        {
            count++;
            continue;
        }

        if (!options.reduce || GlobalV::DRANK == 0)
        {
            detail::write_R_record(ofs,
                                   R_coor.x,
                                   R_coor.y,
                                   R_coor.z,
                                   nonzero_num[count],
                                   options.binary,
                                   text_csr);
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
