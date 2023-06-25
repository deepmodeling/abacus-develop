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
    const double &sparse_threshold, 
    const std::map<Abfs::Vector3_Order<int>,
    std::map<size_t, std::map<size_t, double>>> &pvpR_sparseMatrix,
    LCAO_Matrix *LM
)
{
    ModuleBase::TITLE("Gint_k","distribute_pvpR_sparseMatrix");

    int total_R_num = LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    int *minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : LM->all_R_coor)
    {
        auto iter = pvpR_sparseMatrix.find(R_coor);
        if (iter != pvpR_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvpR_sparseMatrix.find(minus_R_coor);
        if (iter != pvpR_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_int_all(minus_nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(minus_nonzero_num, total_R_num);

    double* tmp = nullptr;
    tmp = new double[GlobalV::NLOCAL];

    count = 0;
    for (auto &R_coor : LM->all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for(int row = 0; row < GlobalV::NLOCAL; ++row)
            {        
                ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
                
                auto iter = pvpR_sparseMatrix.find(R_coor);
                if (iter != pvpR_sparseMatrix.end())
                {
                    
                    if(this->gridt->trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto &value : row_iter->second)
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
                        if(this->gridt->trace_lo[col] >= 0)
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
                
                Parallel_Reduce::reduce_double_pool(tmp, GlobalV::NLOCAL);

                if (LM->ParaV->trace_loc_row[row] >= 0)
                {
                    for(int col = 0; col < GlobalV::NLOCAL; ++col)
                    {
                        if(LM->ParaV->trace_loc_col[col] >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_threshold)
                            {
                                double &value = LM->HR_sparse[current_spin][R_coor][row][col];
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

void Gint_k::distribute_pvpR_soc_sparseMatrix(
    const double &sparse_threshold, 
    const std::map<Abfs::Vector3_Order<int>,
    std::map<size_t, std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix,
    LCAO_Matrix *LM
)
{
    ModuleBase::TITLE("Gint_k","distribute_pvpR_soc_sparseMatrix");

    int total_R_num = LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    int *minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : LM->all_R_coor)
    {
        auto iter = pvpR_soc_sparseMatrix.find(R_coor);
        if (iter != pvpR_soc_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvpR_soc_sparseMatrix.find(minus_R_coor);
        if (iter != pvpR_soc_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_int_all(minus_nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(minus_nonzero_num, total_R_num);

    std::complex<double>* tmp_soc = nullptr;
    tmp_soc = new std::complex<double>[GlobalV::NLOCAL];

    count = 0;
    for (auto &R_coor : LM->all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for(int row = 0; row < GlobalV::NLOCAL; ++row)
            {        
                ModuleBase::GlobalFunc::ZEROS(tmp_soc, GlobalV::NLOCAL);
                
                auto iter = pvpR_soc_sparseMatrix.find(R_coor);
                if (iter != pvpR_soc_sparseMatrix.end())
                {
                    if(this->gridt->trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto &value : row_iter->second)
                            {
                                tmp_soc[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvpR_soc_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvpR_soc_sparseMatrix.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if(this->gridt->trace_lo[col] >= 0)
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
                
                Parallel_Reduce::reduce_complex_double_pool(tmp_soc, GlobalV::NLOCAL);

                if (LM->ParaV->trace_loc_row[row] >= 0)
                {
                    for(int col = 0; col < GlobalV::NLOCAL; ++col)
                    {
                        if(LM->ParaV->trace_loc_col[col] >= 0)
                        {
                            if (std::abs(tmp_soc[col]) > sparse_threshold)
                            {
                                std::complex<double> &value = LM->HR_soc_sparse[R_coor][row][col];
                                value += tmp_soc[col];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    LM->HR_soc_sparse[R_coor][row].erase(col);
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

void Gint_k::cal_vlocal_R_sparseMatrix(const int &current_spin, const double &sparse_threshold, LCAO_Matrix *LM)
{
    ModuleBase::TITLE("Gint_k","cal_vlocal_R_sparseMatrix");

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvpR_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> pvpR_soc_sparseMatrix;

    int lgd = 0;
    double temp_value_double;
    std::complex<double> temp_value_complex;

    ModuleBase::Vector3<double> tau1, dtau;
    for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        for(int I1=0; I1<GlobalC::ucell.atoms[T1].na; ++I1)
        {
            const int iat = GlobalC::ucell.itia2iat(T1,I1);
            if(this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

                const int DM_start = this->gridt->nlocstartg[iat];
                tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
                int nad2 = 0;

                for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
                {
                    const int T2 = GlobalC::GridD.getType(ad);
                    const int I2 = GlobalC::GridD.getNatom(ad);
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);

                    if(this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                        if(distance < rcut)
                        {
                            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                            Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);
                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad2];
                            for(int iw=0; iw<atom1->nw * GlobalV::NPOL; iw++)
                            {
                                for(int iw2=0;iw2<atom2->nw * GlobalV::NPOL; iw2++)
                                {
                                    const int nw = atom2->nw;
                                    const int mug0 = iw/GlobalV::NPOL;
                                    const int nug0 = iw2/GlobalV::NPOL;
                                    const int iw_nowg = ixxx + mug0*nw + nug0;

                                    if(GlobalV::NSPIN == 4)
                                    {		
                                        // pvp is symmetric, only half is calculated.

                                        if(iw%2==0&&iw2%2==0)
                                        {
                                            //spin = 0;
                                            temp_value_complex = std::complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] + std::complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
                                            if(std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                            }
                                        }	
                                        else if(iw%2==1&&iw2%2==1)
                                        {
                                            //spin = 3;
                                            temp_value_complex = std::complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] - std::complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
                                            if(std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                            }
                                        }
                                        else if(iw%2==0&&iw2%2==1)
                                        {
                                            // spin = 1;
                                            if(!GlobalV::DOMAG)
                                            {
                                                // do nothing
                                            }
                                            else
                                            {
                                                temp_value_complex = pvpR_reduced[1][iw_nowg] - std::complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
                                                if(std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                                }
                                            }
                                        }	
                                        else if(iw%2==1&&iw2%2==0) 
                                        {
                                            //spin = 2;
                                            if(!GlobalV::DOMAG)
                                            {
                                                // do nothing
                                            }
                                            else
                                            {
                                                temp_value_complex = pvpR_reduced[1][iw_nowg] + std::complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
                                                if(std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                                }
                                            }
                                        }
                                        else
                                        {
                                            ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                                        }
                                    } //endif NC
                                    else
                                    {
                                        temp_value_double = pvpR_reduced[current_spin][iw_nowg];
                                        if (std::abs(temp_value_double) > sparse_threshold)
                                        {
                                            pvpR_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_double;
                                        }

                                    } //endif normal

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

    if (GlobalV::NSPIN != 4)
    {
        distribute_pvpR_sparseMatrix(current_spin, sparse_threshold, pvpR_sparseMatrix, LM);
    }
    else
    {
        distribute_pvpR_soc_sparseMatrix(sparse_threshold, pvpR_soc_sparseMatrix, LM);
    }

    return;
}