#include "gint_k.h"
#include "grid_technique.h"
#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_k::distribute_pvdpR_sparseMatrix(
    const int current_spin,
    const int dim,
    const double& sparse_threshold,
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& pvdpR_sparseMatrix,
    LCAO_HS_Arrays& HS_Arrays,
    const Parallel_Orbitals* pv)
{
    ModuleBase::TITLE("Gint_k", "distribute_pvdpR_sparseMatrix");

    int total_R_num = HS_Arrays.all_R_coor.size();
    int* nonzero_num = new int[total_R_num];
    int* minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto& R_coor: HS_Arrays.all_R_coor)
    {
        auto iter = pvdpR_sparseMatrix.find(R_coor);
        if (iter != pvdpR_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvdpR_sparseMatrix.find(minus_R_coor);
        if (iter != pvdpR_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_all(minus_nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_pool(nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_pool(minus_nonzero_num, total_R_num);

    double* tmp = nullptr;
    tmp = new double[PARAM.globalv.nlocal];

    count = 0;
    for (auto& R_coor: HS_Arrays.all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for (int row = 0; row < PARAM.globalv.nlocal; ++row)
            {
                ModuleBase::GlobalFunc::ZEROS(tmp, PARAM.globalv.nlocal);

                auto iter = pvdpR_sparseMatrix.find(R_coor);
                if (iter != pvdpR_sparseMatrix.end())
                {

                    if (this->gridt->trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto& value: row_iter->second)
                            {
                                tmp[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvdpR_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvdpR_sparseMatrix.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if (this->gridt->trace_lo[col] >= 0)
                        {
                            auto row_iter = minus_R_iter->second.find(col);
                            if (row_iter != minus_R_iter->second.end())
                            {
                                auto col_iter = row_iter->second.find(row);
                                if (col_iter != row_iter->second.end())
                                {
                                    tmp[col] = col_iter->second;
                                }
                            }
                        }
                    }
                }

                Parallel_Reduce::reduce_pool(tmp, PARAM.globalv.nlocal);

                if (pv->global2local_row(row) >= 0)
                {
                    for (int col = 0; col < PARAM.globalv.nlocal; ++col)
                    {
                        if (pv->global2local_col(col) >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_threshold)
                            {
                                if (dim == 0)
                                {
                                    double& value = HS_Arrays.dHRx_sparse[current_spin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        HS_Arrays.dHRx_sparse[current_spin][R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 1)
                                {
                                    double& value = HS_Arrays.dHRy_sparse[current_spin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        HS_Arrays.dHRy_sparse[current_spin][R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 2)
                                {
                                    double& value = HS_Arrays.dHRz_sparse[current_spin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        HS_Arrays.dHRz_sparse[current_spin][R_coor][row].erase(col);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        count++;
    }

    delete[] nonzero_num;
    delete[] minus_nonzero_num;
    delete[] tmp;
    nonzero_num = nullptr;
    minus_nonzero_num = nullptr;
    tmp = nullptr;

    return;
}

void Gint_k::distribute_pvdpR_soc_sparseMatrix(
    const int dim,
    const double& sparse_threshold,
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>&
        pvdpR_soc_sparseMatrix,
    LCAO_HS_Arrays& HS_Arrays,
    const Parallel_Orbitals* pv)
{
    ModuleBase::TITLE("Gint_k", "distribute_pvdpR_soc_sparseMatrix");

    int total_R_num = HS_Arrays.all_R_coor.size();
    int* nonzero_num = new int[total_R_num];
    int* minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto& R_coor: HS_Arrays.all_R_coor)
    {
        auto iter = pvdpR_soc_sparseMatrix.find(R_coor);
        if (iter != pvdpR_soc_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvdpR_soc_sparseMatrix.find(minus_R_coor);
        if (iter != pvdpR_soc_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_all(minus_nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_pool(nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_pool(minus_nonzero_num, total_R_num);

    std::complex<double>* tmp_soc = nullptr;
    tmp_soc = new std::complex<double>[PARAM.globalv.nlocal];

    count = 0;
    for (auto& R_coor: HS_Arrays.all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for (int row = 0; row < PARAM.globalv.nlocal; ++row)
            {
                ModuleBase::GlobalFunc::ZEROS(tmp_soc, PARAM.globalv.nlocal);

                auto iter = pvdpR_soc_sparseMatrix.find(R_coor);
                if (iter != pvdpR_soc_sparseMatrix.end())
                {
                    if (this->gridt->trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto& value: row_iter->second)
                            {
                                tmp_soc[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvdpR_soc_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvdpR_soc_sparseMatrix.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if (this->gridt->trace_lo[col] >= 0)
                        {
                            auto row_iter = minus_R_iter->second.find(col);
                            if (row_iter != minus_R_iter->second.end())
                            {
                                auto col_iter = row_iter->second.find(row);
                                if (col_iter != row_iter->second.end())
                                {
                                    tmp_soc[col] = conj(col_iter->second);
                                }
                            }
                        }
                    }
                }

                Parallel_Reduce::reduce_pool(tmp_soc, PARAM.globalv.nlocal);

                if (pv->global2local_row(row) >= 0)
                {
                    for (int col = 0; col < PARAM.globalv.nlocal; ++col)
                    {
                        if (pv->global2local_col(col) >= 0)
                        {
                            if (std::abs(tmp_soc[col]) > sparse_threshold)
                            {
                                if (dim == 0)
                                {
                                    std::complex<double>& value = HS_Arrays.dHRx_soc_sparse[R_coor][row][col];
                                    value += tmp_soc[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        HS_Arrays.dHRx_soc_sparse[R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 1)
                                {
                                    std::complex<double>& value = HS_Arrays.dHRy_soc_sparse[R_coor][row][col];
                                    value += tmp_soc[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        HS_Arrays.dHRy_soc_sparse[R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 2)
                                {
                                    std::complex<double>& value = HS_Arrays.dHRz_soc_sparse[R_coor][row][col];
                                    value += tmp_soc[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        HS_Arrays.dHRz_soc_sparse[R_coor][row].erase(col);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        count++;
    }

    delete[] nonzero_num;
    delete[] minus_nonzero_num;
    delete[] tmp_soc;
    nonzero_num = nullptr;
    minus_nonzero_num = nullptr;
    tmp_soc = nullptr;

    return;
}

void Gint_k::cal_dvlocal_R_sparseMatrix(const int& current_spin,
                                        const double& sparse_threshold,
                                        LCAO_HS_Arrays& HS_Arrays,
                                        const Parallel_Orbitals* pv,
                                        UnitCell& ucell,
                                        Grid_Driver& gdriver)
{
    ModuleBase::TITLE("Gint_k", "cal_dvlocal_R_sparseMatrix");

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvdpRx_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvdpRy_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvdpRz_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        pvdpRx_soc_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        pvdpRy_soc_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        pvdpRz_soc_sparseMatrix;

    int lgd = 0;
    double temp_value_double;
    std::complex<double> temp_value_complex;

    ModuleBase::Vector3<double> tau1, dtau;
    for (int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        for (int I1 = 0; I1 < ucell.atoms[T1].na; ++I1)
        {
            const int iat = ucell.itia2iat(T1, I1);
            if (this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &ucell.atoms[T1];
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

                const int DM_start = this->gridt->nlocstartg[iat];
                tau1 = ucell.atoms[T1].tau[I1];
                gdriver.Find_atom(ucell, tau1, T1, I1);
                int nad2 = 0;

                for (int ad = 0; ad < gdriver.getAdjacentNum() + 1; ad++)
                {
                    const int T2 = gdriver.getType(ad);
                    const int I2 = gdriver.getNatom(ad);
                    const int iat2 = ucell.itia2iat(T2, I2);

                    if (this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &ucell.atoms[T2];
                        dtau = gdriver.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * ucell.lat0;
                        double rcut = this->gridt->rcuts[T1] + this->gridt->rcuts[T2];

                        if (distance < rcut)
                        {
                            const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                            Abfs::Vector3_Order<int> dR(gdriver.getBox(ad).x,
                                                        gdriver.getBox(ad).y,
                                                        gdriver.getBox(ad).z);

                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad2];

                            for (int iw = 0; iw < atom1->nw * PARAM.globalv.npol; iw++)
                            {
                                for (int iw2 = 0; iw2 < atom2->nw * PARAM.globalv.npol; iw2++)
                                {
                                    const int nw = atom2->nw;
                                    const int mug0 = iw / PARAM.globalv.npol;
                                    const int nug0 = iw2 / PARAM.globalv.npol;
                                    const int iw_nowg = ixxx + mug0 * nw + nug0;

                                    if (PARAM.inp.nspin == 4)
                                    {
                                        // pvp is symmetric, only half is calculated.

                                        if (iw % 2 == 0 && iw2 % 2 == 0)
                                        {
                                            // spin = 0;
                                            temp_value_complex
                                                = std::complex<double>(1.0, 0.0) * pvdpRx_reduced[0][iw_nowg]
                                                  + std::complex<double>(1.0, 0.0) * pvdpRx_reduced[3][iw_nowg];

                                            if (std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvdpRx_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                    = temp_value_complex;
                                            }

                                            temp_value_complex
                                                = std::complex<double>(1.0, 0.0) * pvdpRy_reduced[0][iw_nowg]
                                                  + std::complex<double>(1.0, 0.0) * pvdpRy_reduced[3][iw_nowg];

                                            if (std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvdpRy_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                    = temp_value_complex;
                                            }
                                            temp_value_complex
                                                = std::complex<double>(1.0, 0.0) * pvdpRz_reduced[0][iw_nowg]
                                                  + std::complex<double>(1.0, 0.0) * pvdpRz_reduced[3][iw_nowg];

                                            if (std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvdpRz_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                    = temp_value_complex;
                                            }
                                        }
                                        else if (iw % 2 == 1 && iw2 % 2 == 1)
                                        {
                                            // spin = 3;
                                            temp_value_complex
                                                = std::complex<double>(1.0, 0.0) * pvdpRx_reduced[0][iw_nowg]
                                                  - std::complex<double>(1.0, 0.0) * pvdpRx_reduced[3][iw_nowg];
                                            if (std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvdpRx_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                    = temp_value_complex;
                                            }
                                            temp_value_complex
                                                = std::complex<double>(1.0, 0.0) * pvdpRy_reduced[0][iw_nowg]
                                                  - std::complex<double>(1.0, 0.0) * pvdpRy_reduced[3][iw_nowg];
                                            if (std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvdpRy_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                    = temp_value_complex;
                                            }
                                            temp_value_complex
                                                = std::complex<double>(1.0, 0.0) * pvdpRz_reduced[0][iw_nowg]
                                                  - std::complex<double>(1.0, 0.0) * pvdpRz_reduced[3][iw_nowg];
                                            if (std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvdpRz_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                    = temp_value_complex;
                                            }
                                        }
                                        else if (iw % 2 == 0 && iw2 % 2 == 1)
                                        {
                                            // spin = 1;
                                            if (!PARAM.globalv.domag)
                                            {
                                                // do nothing
                                            }
                                            else
                                            {
                                                temp_value_complex
                                                    = pvdpRx_reduced[1][iw_nowg]
                                                      - std::complex<double>(0.0, 1.0) * pvdpRx_reduced[2][iw_nowg];
                                                if (std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvdpRx_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                        = temp_value_complex;
                                                }
                                                temp_value_complex
                                                    = pvdpRy_reduced[1][iw_nowg]
                                                      - std::complex<double>(0.0, 1.0) * pvdpRy_reduced[2][iw_nowg];
                                                if (std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvdpRy_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                        = temp_value_complex;
                                                }
                                                temp_value_complex
                                                    = pvdpRz_reduced[1][iw_nowg]
                                                      - std::complex<double>(0.0, 1.0) * pvdpRz_reduced[2][iw_nowg];
                                                if (std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvdpRz_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                        = temp_value_complex;
                                                }
                                            }
                                        }
                                        else if (iw % 2 == 1 && iw2 % 2 == 0)
                                        {
                                            // spin = 2;
                                            if (!PARAM.globalv.domag)
                                            {
                                                // do nothing
                                            }
                                            else
                                            {
                                                temp_value_complex
                                                    = pvdpRx_reduced[1][iw_nowg]
                                                      + std::complex<double>(0.0, 1.0) * pvdpRx_reduced[2][iw_nowg];
                                                if (std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvdpRx_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                        = temp_value_complex;
                                                }
                                                temp_value_complex
                                                    = pvdpRy_reduced[1][iw_nowg]
                                                      + std::complex<double>(0.0, 1.0) * pvdpRy_reduced[2][iw_nowg];
                                                if (std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvdpRy_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                        = temp_value_complex;
                                                }
                                                temp_value_complex
                                                    = pvdpRz_reduced[1][iw_nowg]
                                                      + std::complex<double>(0.0, 1.0) * pvdpRz_reduced[2][iw_nowg];
                                                if (std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvdpRz_soc_sparseMatrix[dR][start1 + iw][start2 + iw2]
                                                        = temp_value_complex;
                                                }
                                            }
                                        }
                                        else
                                        {
                                            ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc", "index is wrong!");
                                        }
                                    } // endif NC
                                    else
                                    {
                                        temp_value_double = pvdpRx_reduced[current_spin][iw_nowg];
                                        if (std::abs(temp_value_double) > sparse_threshold)
                                        {
                                            pvdpRx_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_double;
                                        }
                                        temp_value_double = pvdpRy_reduced[current_spin][iw_nowg];
                                        if (std::abs(temp_value_double) > sparse_threshold)
                                        {
                                            pvdpRy_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_double;
                                        }
                                        temp_value_double = pvdpRz_reduced[current_spin][iw_nowg];
                                        if (std::abs(temp_value_double) > sparse_threshold)
                                        {
                                            pvdpRz_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_double;
                                        }
                                    } // endif normal
                                }

                                ++lgd;
                            }
                            ++nad2;
                        }
                    }
                }
            }
        }
    }

    if (PARAM.inp.nspin != 4)
    {
        distribute_pvdpR_sparseMatrix(current_spin, 0, sparse_threshold, pvdpRx_sparseMatrix, HS_Arrays, pv);
        distribute_pvdpR_sparseMatrix(current_spin, 1, sparse_threshold, pvdpRy_sparseMatrix, HS_Arrays, pv);
        distribute_pvdpR_sparseMatrix(current_spin, 2, sparse_threshold, pvdpRz_sparseMatrix, HS_Arrays, pv);
    }
    else
    {
        distribute_pvdpR_soc_sparseMatrix(0, sparse_threshold, pvdpRx_soc_sparseMatrix, HS_Arrays, pv);
        distribute_pvdpR_soc_sparseMatrix(1, sparse_threshold, pvdpRy_soc_sparseMatrix, HS_Arrays, pv);
        distribute_pvdpR_soc_sparseMatrix(2, sparse_threshold, pvdpRz_soc_sparseMatrix, HS_Arrays, pv);
    }

    return;
}
