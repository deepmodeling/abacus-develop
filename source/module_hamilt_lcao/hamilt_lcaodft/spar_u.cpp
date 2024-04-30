#include "spar_u.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"

void sparse_format::cal_HR_dftu(
	    const Parallel_Orbitals &pv,
        std::set<Abfs::Vector3_Order<int>> &all_R_coor,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &SR_sparse,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> *HR_sparse,
		const int &current_spin, 
		const double &sparse_thr)
{
    ModuleBase::TITLE("sparse_format","cal_HR_dftu");
    ModuleBase::timer::tick("sparse_format","cal_HR_dftu");

    int total_R_num = all_R_coor.size();
    int *nonzero_num = new int[total_R_num]();
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);

    int count = 0;
    for (auto &R_coor : all_R_coor)
    {
        auto iter = SR_sparse.find(R_coor);
        if (iter != SR_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);

    double *HR_tmp = new double[pv.nloc];
    double *SR_tmp = new double[pv.nloc];

    int ir=0;
    int ic=0;
    int iic=0;
    auto &temp_HR_sparse = HR_sparse[current_spin];

    count = 0;
    for (auto &R_coor : all_R_coor)
    {
        if (nonzero_num[count] != 0)
        {
            ModuleBase::GlobalFunc::ZEROS(HR_tmp, pv.nloc);
            ModuleBase::GlobalFunc::ZEROS(SR_tmp, pv.nloc);

            auto iter = SR_sparse.find(R_coor);
            if (iter != SR_sparse.end())
            {
                for (auto &row_loop : iter->second)
                {
                    ir = pv.global2local_row(row_loop.first);
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = pv.global2local_col(col_loop.first);
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        SR_tmp[iic] = col_loop.second;
                    }
                }
            }

            GlobalC::dftu.cal_eff_pot_mat_R_double(current_spin, SR_tmp, HR_tmp);

            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                ir = pv.global2local_row(i);
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = pv.global2local_col(j);
                        if (ic >= 0)
                        {
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic = ir + ic * pv.nrow;
                            }
                            else
                            {
                                iic = ir * pv.ncol + ic;
                            }

                            if (std::abs(HR_tmp[iic]) > sparse_thr)
                            {
                                double &value = temp_HR_sparse[R_coor][i][j];
                                value += HR_tmp[iic];
                                if (std::abs(value) <= sparse_thr)
                                {
                                    temp_HR_sparse[R_coor][i].erase(j);
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
    delete[] HR_tmp;
    delete[] SR_tmp;
    nonzero_num = nullptr;
    HR_tmp = nullptr;
    SR_tmp = nullptr;

    ModuleBase::timer::tick("sparse_format","cal_HR_dftu_sparse");

    return;
}


void sparse_format::cal_HR_dftu_soc(
	    const Parallel_Orbitals &pv,
        std::set<Abfs::Vector3_Order<int>> &all_R_coor,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &SR_soc_sparse,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &HR_soc_sparse,
		const int &current_spin, 
		const double &sparse_thr)
{
    ModuleBase::TITLE("sparse_format","cal_HR_dftu_soc");
    ModuleBase::timer::tick("sparse_format","cal_HR_dftu_soc");

    int total_R_num = all_R_coor.size();
    int *nonzero_num = new int[total_R_num]();
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : all_R_coor)
    {
        auto iter = SR_soc_sparse.find(R_coor);
        if (iter != SR_soc_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);

    std::complex<double> *HR_soc_tmp = new std::complex<double>[pv.nloc];
    std::complex<double> *SR_soc_tmp = new std::complex<double>[pv.nloc];

    int ir=0;
    int ic=0;
    int iic=0;

    count = 0;
    for (auto &R_coor : all_R_coor)
    {
        if (nonzero_num[count] != 0)
        {
            ModuleBase::GlobalFunc::ZEROS(HR_soc_tmp, pv.nloc);
            ModuleBase::GlobalFunc::ZEROS(SR_soc_tmp, pv.nloc);

            auto iter = SR_soc_sparse.find(R_coor);
            if (iter != SR_soc_sparse.end())
            {
                for (auto &row_loop : iter->second)
                {
                    ir = pv.global2local_row(row_loop.first);
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = pv.global2local_col(col_loop.first);
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        SR_soc_tmp[iic] = col_loop.second;
                    }
                }
            }

            GlobalC::dftu.cal_eff_pot_mat_R_complex_double(current_spin, SR_soc_tmp, HR_soc_tmp);

            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                ir = pv.global2local_row(i);
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = pv.global2local_col(j);
                        if (ic >= 0)
                        {
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic = ir + ic * pv.nrow;
                            }
                            else
                            {
                                iic = ir * pv.ncol + ic;
                            }

                            if (std::abs(HR_soc_tmp[iic]) > sparse_thr)
                            {
                                std::complex<double> &value = HR_soc_sparse[R_coor][i][j];
                                value += HR_soc_tmp[iic];
                                if (std::abs(value) <= sparse_thr)
                                {
                                    HR_soc_sparse[R_coor][i].erase(j);
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
    delete[] HR_soc_tmp;
    delete[] SR_soc_tmp;
    nonzero_num = nullptr;
    HR_soc_tmp = nullptr;
    SR_soc_tmp = nullptr;

    ModuleBase::timer::tick("sparse_format","calculat_HR_dftu_soc");

    return;
}
