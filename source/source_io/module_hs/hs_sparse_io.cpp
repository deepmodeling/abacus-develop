#include "hs_sparse_io.h"
#include "hs_sparse_io_detail.h"

#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
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
    const bool append_on_restart = (options.calculation == "md") && options.out_app_flag && options.istep >= 0;
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
    const int step = std::max(options.istep, 0);
    if (options.binary)
    {
        ofs.write(reinterpret_cast<const char*>(&step), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&nlocal), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&output_R_number), sizeof(int));
    }
    else
    {
        ofs << "STEP: " << step << std::endl;
        ofs << "Matrix Dimension of " + options.label + "(R): " << nlocal
            << std::endl;
        ofs << "Matrix number of " + options.label + "(R): "
            << output_R_number << std::endl;
    }
}

void write_R_record(std::ofstream& ofs,
                    const ModuleIO::RCoordinate& R_coor,
                    const long long nonzero_count,
                    const bool binary)
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
        ofs << dRx << " " << dRy << " " << dRz << " " << nonzero_count
            << std::endl;
    }
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
    const int nlocal = pv.get_global_row_size();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_sparse",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    const std::vector<long long> nonzero_num
        = detail::count_nonzeros_by_R(smat, all_R_coor, options.threshold, options.reduce);
    const int output_R_number = count_output_R(nonzero_num);
    std::ofstream ofs;
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        open_sparse_file(ofs, options);
        write_sparse_header(ofs, options, nlocal, output_R_number);
    }

    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        if (nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (!options.reduce || GlobalV::DRANK == 0)
        {
            write_R_record(ofs, R_coor, nonzero_num[count], options.binary);
        }

        if (smat.count(R_coor))
        {
            save_lat_r(ofs, smat.at(R_coor), pv, options);
        }
        else
        {
            SparseRBlock<Tdata> empty_map;
            save_lat_r(ofs, empty_map, pv, options);
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
