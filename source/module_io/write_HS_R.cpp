#include "write_HS_R.h"

#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_dh.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_hsr.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_st.h"
#include "write_HS_sparse.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <vector>

namespace
{
template <typename T>
struct GetSEntry
{
    size_t row;
    size_t col;
    T value;
};

void write_gets_data(std::ofstream& ofs, const double& data)
{
    ofs << " " << std::fixed << std::scientific << std::setprecision(8) << data;
}

void write_gets_data(std::ofstream& ofs, const std::complex<double>& data)
{
    ofs << " (" << std::fixed << std::scientific << std::setprecision(8) << data.real() << ","
        << std::fixed << std::scientific << std::setprecision(8) << data.imag() << ")";
}

template <typename Tdata>
void save_gets_sparse_fast(
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>>& smat,
    const std::set<Abfs::Vector3_Order<int>>& all_R_coor,
    const double& sparse_thr,
    const std::string& filename,
    const std::string& label,
    const int& istep)
{
    const int nlocal = PARAM.globalv.nlocal;
    const int total_R_num = all_R_coor.size();

    std::vector<int> nonzero_num(total_R_num, 0);
    int count = 0;
    for (const auto& R_coor : all_R_coor)
    {
        auto iter = smat.find(R_coor);
        if (iter != smat.end())
        {
            for (const auto& row_loop : iter->second)
            {
                for (const auto& value : row_loop.second)
                {
                    if (std::abs(value.second) > sparse_thr)
                    {
                        ++nonzero_num[count];
                    }
                }
            }
        }
        ++count;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(nonzero_num.data(), total_R_num);
#endif

    int output_R_number = 0;
    for (const int nnz : nonzero_num)
    {
        if (nnz != 0)
        {
            ++output_R_number;
        }
    }

    int myrank = 0;
    int nprocs = 1;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

    std::ofstream ofs;
    if (myrank == 0)
    {
        ofs.open(filename.c_str());
        ofs << "STEP: " << std::max(istep, 0) << std::endl;
        ofs << "Matrix Dimension of " + label + "(R): " << nlocal << std::endl;
        ofs << "Matrix number of " + label + "(R): " << output_R_number << std::endl;
    }

    count = 0;
    for (const auto& R_coor : all_R_coor)
    {
        if (nonzero_num[count] == 0)
        {
            ++count;
            continue;
        }

        std::ofstream ofs_indices;
        std::ifstream ifs_indices;
        std::stringstream indices_filename;
        indices_filename << PARAM.globalv.global_out_dir << std::to_string(myrank) << "temp_gets_sparse_indices.dat";

        std::vector<int> indptr;
        if (myrank == 0)
        {
            ofs << R_coor.x << " " << R_coor.y << " " << R_coor.z << " "
                << nonzero_num[count] << std::endl;
            ofs_indices.open(indices_filename.str().c_str());
            indptr.reserve(nlocal + 1);
            indptr.push_back(0);
        }

        const int row_block_size = 4096;
        for (int row_begin = 0; row_begin < nlocal; row_begin += row_block_size)
        {
            const int row_end = std::min(row_begin + row_block_size, nlocal);

            std::vector<GetSEntry<Tdata>> local_entries;
            auto iter = smat.find(R_coor);
            if (iter != smat.end())
            {
                auto row_iter = iter->second.lower_bound(static_cast<size_t>(row_begin));
                while (row_iter != iter->second.end() && row_iter->first < static_cast<size_t>(row_end))
                {
                    for (const auto& value : row_iter->second)
                    {
                        if (std::abs(value.second) > sparse_thr)
                        {
                            local_entries.push_back({row_iter->first, value.first, value.second});
                        }
                    }
                    ++row_iter;
                }
            }

            std::vector<GetSEntry<Tdata>> entries;
#ifdef __MPI
            const int local_bytes = static_cast<int>(local_entries.size() * sizeof(GetSEntry<Tdata>));
            std::vector<int> recv_counts(nprocs, 0);
            MPI_Gather(&local_bytes, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

            std::vector<int> displs;
            int total_bytes = 0;
            if (myrank == 0)
            {
                displs.resize(nprocs, 0);
                for (int ip = 1; ip < nprocs; ++ip)
                {
                    displs[ip] = displs[ip - 1] + recv_counts[ip - 1];
                }
                total_bytes = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
                entries.resize(total_bytes / sizeof(GetSEntry<Tdata>));
            }

            MPI_Gatherv(local_entries.empty() ? nullptr : local_entries.data(),
                        local_bytes,
                        MPI_BYTE,
                        entries.empty() ? nullptr : entries.data(),
                        recv_counts.empty() ? nullptr : recv_counts.data(),
                        displs.empty() ? nullptr : displs.data(),
                        MPI_BYTE,
                        0,
                        MPI_COMM_WORLD);
#else
            entries.swap(local_entries);
#endif

            if (myrank == 0)
            {
                std::sort(entries.begin(), entries.end(), [](const auto& lhs, const auto& rhs) {
                    return lhs.row == rhs.row ? lhs.col < rhs.col : lhs.row < rhs.row;
                });

                std::vector<GetSEntry<Tdata>> merged_entries;
                merged_entries.reserve(entries.size());
                for (size_t i = 0; i < entries.size();)
                {
                    GetSEntry<Tdata> merged = entries[i];
                    ++i;
                    while (i < entries.size() && entries[i].row == merged.row && entries[i].col == merged.col)
                    {
                        merged.value += entries[i].value;
                        ++i;
                    }
                    if (std::abs(merged.value) > sparse_thr)
                    {
                        merged_entries.push_back(merged);
                    }
                }

                size_t entry_index = 0;
                for (int row = row_begin; row < row_end; ++row)
                {
                    int row_nnz = 0;
                    while (entry_index < merged_entries.size()
                           && merged_entries[entry_index].row == static_cast<size_t>(row))
                    {
                        write_gets_data(ofs, merged_entries[entry_index].value);
                        ofs_indices << " " << merged_entries[entry_index].col;
                        ++row_nnz;
                        ++entry_index;
                    }
                    indptr.push_back(indptr.back() + row_nnz);
                }
            }
        }

        if (myrank == 0)
        {
            ofs << std::endl;

            ofs_indices << std::endl;
            ofs_indices.close();
            ifs_indices.open(indices_filename.str().c_str());
            ofs << ifs_indices.rdbuf();
            ifs_indices.close();

            for (const auto& pointer : indptr)
            {
                ofs << " " << pointer;
            }
            ofs << std::endl;

            std::remove(indices_filename.str().c_str());
        }

        ++count;
    }

    if (myrank == 0)
    {
        ofs.close();
    }
}
} // namespace

// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix.
// If the absolute value of the matrix element is less than or equal to the
// 'sparse_thr', it will be ignored.
void ModuleIO::output_HSR(const UnitCell& ucell,
                          const int& istep,
                          const ModuleBase::matrix& v_eff,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const K_Vectors& kv,
                          hamilt::Hamilt<std::complex<double>>* p_ham,
#ifdef __EXX
                          const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd,
                          const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc,
#endif
                          const std::string& SR_filename,
                          const std::string& HR_filename_up,
                          const std::string HR_filename_down,
                          const bool& binary,
                          const double& sparse_thr) {
    ModuleBase::TITLE("ModuleIO", "output_HSR");
    ModuleBase::timer::tick("ModuleIO", "output_HSR");

    const int nspin = PARAM.inp.nspin;

    if (nspin == 1 || nspin == 4) {
        const int spin_now = 0;
        // jingan add 2021-6-4, modify 2021-12-2
        sparse_format::cal_HSR(ucell,pv, HS_Arrays, grid, spin_now, sparse_thr, kv.nmp, p_ham
#ifdef __EXX
            , Hexxd, Hexxc
#endif
        );
    }
    else if (nspin == 2) {
        int spin_now = 1;

        // save HR of spin down first (the current spin always be down)
        sparse_format::cal_HSR(ucell,pv, HS_Arrays, grid, spin_now, sparse_thr, kv.nmp, p_ham
#ifdef __EXX
            , Hexxd, Hexxc
#endif
        );

        // cal HR of the spin up
        if (PARAM.inp.vl_in_h) {
            const int ik = 0;
            p_ham->refresh();
            p_ham->updateHk(ik);
            spin_now = 0;
        }

        sparse_format::cal_HSR(ucell,pv, HS_Arrays, grid, spin_now, sparse_thr, kv.nmp, p_ham
#ifdef __EXX
            , Hexxd, Hexxc
#endif
        );
    }

    ModuleIO::save_HSR_sparse(istep,
                              pv,
                              HS_Arrays,
                              sparse_thr,
                              binary,
                              SR_filename,
                              HR_filename_up,
                              HR_filename_down);

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_HSR");
    return;
}

void ModuleIO::output_dHR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          Gint_k& gint_k, // mohan add 2024-04-01
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_dHR");
    ModuleBase::timer::tick("ModuleIO", "output_dHR");

    gint_k.allocate_pvdpR();

    const int nspin = PARAM.inp.nspin;

    if (nspin == 1 || nspin == 4) {
        // mohan add 2024-04-01
        const int cspin = 0;

        sparse_format::cal_dH(ucell,
                              pv,
                              HS_Arrays,
                              grid,
                              two_center_bundle,
                              orb,
                              cspin,
                              sparse_thr,
                              gint_k);
    } else if (nspin == 2) {
        for (int cspin = 0; cspin < 2; cspin++) {
            // note: some MPI process will not have grids when MPI cores are too
            // many, v_eff in these processes are empty
            const double* vr_eff1
                = v_eff.nc * v_eff.nr > 0 ? &(v_eff(cspin, 0)) : nullptr;

            if (!PARAM.globalv.gamma_only_local) {
                if (PARAM.inp.vl_in_h) {
                    Gint_inout inout(vr_eff1,
                                     cspin,
                                     Gint_Tools::job_type::dvlocal);
                    gint_k.cal_gint(&inout);
                }
            }

            sparse_format::cal_dH(ucell,
                                  pv,
                                  HS_Arrays,
                                  grid,
                                  two_center_bundle,
                                  orb,
                                  cspin,
                                  sparse_thr,
                                  gint_k);
        }
    }
    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary);

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    gint_k.destroy_pvdpR();

    ModuleBase::timer::tick("ModuleIO", "output_dHR");
    return;
}

void ModuleIO::output_SR(Parallel_Orbitals& pv,
                         const Grid_Driver& grid,
                         hamilt::Hamilt<std::complex<double>>* p_ham,
                         const std::string& SR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_SR");
    ModuleBase::timer::tick("ModuleIO", "output_SR");

    LCAO_HS_Arrays HS_Arrays;

    sparse_format::cal_SR(pv,
                          HS_Arrays.all_R_coor,
                          HS_Arrays.SR_sparse,
                          HS_Arrays.SR_soc_sparse,
                          grid,
                          sparse_thr,
                          p_ham);

    const int istep = 0;

    if (PARAM.inp.nspin == 4)
    {
        if (PARAM.inp.calculation == "get_S" && !binary)
        {
            save_gets_sparse_fast(HS_Arrays.SR_soc_sparse,
                                  HS_Arrays.all_R_coor,
                                  sparse_thr,
                                  SR_filename,
                                  "S",
                                  istep);
        }
        else
        {
            ModuleIO::save_sparse(HS_Arrays.SR_soc_sparse,
                                  HS_Arrays.all_R_coor,
                                  sparse_thr,
                                  binary,
                                  SR_filename,
                                  pv,
                                  "S",
                                  istep);
        }
    }
    else
    {
        if (PARAM.inp.calculation == "get_S" && !binary)
        {
            save_gets_sparse_fast(HS_Arrays.SR_sparse,
                                  HS_Arrays.all_R_coor,
                                  sparse_thr,
                                  SR_filename,
                                  "S",
                                  istep);
        }
        else
        {
            ModuleIO::save_sparse(HS_Arrays.SR_sparse,
                                  HS_Arrays.all_R_coor,
                                  sparse_thr,
                                  binary,
                                  SR_filename,
                                  pv,
                                  "S",
                                  istep);
        }
    }

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_SR");
    return;
}

void ModuleIO::output_TR(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         LCAO_HS_Arrays& HS_Arrays,
                         const Grid_Driver& grid,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const std::string& TR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_TR");
    ModuleBase::timer::tick("ModuleIO", "output_TR");

    std::stringstream sst;
    if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag) {
        sst << PARAM.globalv.global_matrix_dir << istep << "_" << TR_filename;
    } else {
        sst << PARAM.globalv.global_out_dir << TR_filename;
    }

    sparse_format::cal_TR(ucell,
                          pv,
                          HS_Arrays,
                          grid,
                          two_center_bundle,
                          orb,
                          sparse_thr);

    ModuleIO::save_sparse(HS_Arrays.TR_sparse,
                          HS_Arrays.all_R_coor,
                          sparse_thr,
                          binary,
                          sst.str().c_str(),
                          pv,
                          "T",
                          istep);

    sparse_format::destroy_T_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_TR");
    return;
}
