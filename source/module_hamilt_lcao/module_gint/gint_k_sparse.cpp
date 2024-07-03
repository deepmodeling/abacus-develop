#include "gint_k.h"
#include "grid_technique.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_k::distribute_pvpR_sparseMatrix(
    const int current_spin,
    const double& sparse_threshold,
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& pvpR_sparseMatrix,
    LCAO_Matrix* LM,
    Parallel_Orbitals* pv)
{
    ModuleBase::TITLE("Gint_k", "distribute_pvpR_sparseMatrix");

    int total_R_num = LM->all_R_coor.size();
    int* nonzero_num = new int[total_R_num];
    int* minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto& R_coor: LM->all_R_coor)
    {
        auto iter = pvpR_sparseMatrix.find(R_coor);
        if (iter != pvpR_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvpR_sparseMatrix.find(minus_R_coor);
        if (iter != pvpR_sparseMatrix.end())
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
    tmp = new double[GlobalV::NLOCAL];

    count = 0;
    for (auto& R_coor: LM->all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for (int row = 0; row < GlobalV::NLOCAL; ++row)
            {
                ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);

                auto iter = pvpR_sparseMatrix.find(R_coor);
                if (iter != pvpR_sparseMatrix.end())
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

                auto minus_R_iter = pvpR_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvpR_sparseMatrix.end())
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

                Parallel_Reduce::reduce_pool(tmp, GlobalV::NLOCAL);

                if (pv->global2local_row(row) >= 0)
                {
                    for (int col = 0; col < GlobalV::NLOCAL; ++col)
                    {
                        if (pv->global2local_col(col) >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_threshold)
                            {
                                double& value = LM->HR_sparse[current_spin][R_coor][row][col];
                                value += tmp[col];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    LM->HR_sparse[current_spin][R_coor][row].erase(col);
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