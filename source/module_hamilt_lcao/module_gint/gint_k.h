#ifndef GINT_K_H
#define GINT_K_H

#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "grid_technique.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_elecstate/module_charge/charge.h"
#include "gint.h"
#include "module_base/parallel_reduce.h"
#include "module_basis/module_ao/ORB_read.h"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "module_base/abfs-vector3_order.h"
#include <functional>


class Gint_k : public Gint
{
    public:
    ~Gint_k()
    {
        destroy_pvpR();
    }
    //------------------------------------------------------
    // in gint_k_pvpr.cpp 
    //------------------------------------------------------
    // pvpR and reset_spin/get_spin : auxilliary methods
    // for calculating hamiltonian

    // reset the spin.
    void reset_spin(const int &spin_now_in){this->spin_now = spin_now_in;};
    // get the spin.
    int get_spin(void)const{return spin_now;}

    //renew gint index for new iteration
    void renew(const bool& soft = false)
    {
        if(soft && this->spin_now == 0) 
        {//in this case, gint will not be recalculated
            return;
        }
        else if (this->spin_now != -1)
        {
            int start_spin = -1;
            this->reset_spin(start_spin);
            this->destroy_pvpR();
            this->allocate_pvpR();
        }
        return;
    }
 
    // allocate the <phi_0 | V | phi_R> matrix element.
    void allocate_pvpR(void);
    // destroy the temporary <phi_0 | V | phi_R> matrix element.
    void destroy_pvpR(void);

    // allocate the <phi_0 | V | dphi_R> matrix element.
    void allocate_pvdpR(void);
    // destroy the temporary <phi_0 | V | dphi_R> matrix element.
    void destroy_pvdpR(void);

    // folding the < phi_0 | V | phi_R> matrix to 
    // <phi_0i | V | phi_0j>
    // V is (Vl + Vh + Vxc) if no Vna is used,
    // and is (Vna + delta_Vh + Vxc) if Vna is used.
    void folding_vl_k(const int &ik, LCAO_Matrix* LM, Parallel_Orbitals *pv,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                    const UnitCell& ucell,const LCAO_Orbitals& orb,Grid_Driver& gd);

    /**
     * @brief transfer pvpR to this->hRGint
     * then pass this->hRGint to Veff<OperatorLCAO>::hR
    */
    void transfer_pvpR(hamilt::HContainer<double> *hR,const UnitCell* ucell_in,const LCAO_Orbitals& orb,Grid_Driver* gd);
    void transfer_pvpR(hamilt::HContainer<std::complex<double>> *hR,const UnitCell* ucell_in,const LCAO_Orbitals& orb,Grid_Driver* gd);

    //------------------------------------------------------
    // in gint_k_env.cpp 
    //------------------------------------------------------
    // calculate the envelop function via grid integrals
    void cal_env_k(int ik,
                   const std::complex<double>* psi_k,
                   double* rho,
                   const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                   const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                   LCAO_Orbitals &orb,UnitCell &ucell);

    //------------------------------------------------------
    // in gint_k_sparse.cpp 
    //------------------------------------------------------    
    // related to sparse matrix
    // jingan add 2021-6-4, modify 2021-12-2
    void distribute_pvpR_sparseMatrix(
        const int current_spin, 
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t, std::map<size_t, double>>> &pvpR_sparseMatrix,
        LCAO_Matrix *LM,
        Parallel_Orbitals *pv);

    void distribute_pvpR_soc_sparseMatrix(
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t,
        std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix,
        LCAO_Matrix *LM,
        Parallel_Orbitals *pv
        );

    void cal_vlocal_R_sparseMatrix(
        const int &current_spin,
        const double &sparse_threshold,
        LCAO_Matrix *LM,
        Parallel_Orbitals *pv,
        LCAO_Orbitals &orb,UnitCell &ucell,
        Grid_Driver &gdriver);

    //------------------------------------------------------
    // in gint_k_sparse1.cpp 
    //------------------------------------------------------  
    // similar to the above 3, just for the derivative
    void distribute_pvdpR_sparseMatrix(
        const int current_spin, 
        const int dim,
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t, std::map<size_t, double>>> &pvdpR_sparseMatrix,
        LCAO_Matrix *LM,
        Parallel_Orbitals *pv);

    void distribute_pvdpR_soc_sparseMatrix(
        const int dim,
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t, std::map<size_t, std::complex<double>>>> &pvdpR_soc_sparseMatrix,
        LCAO_Matrix *LM,
        Parallel_Orbitals *pv);

    void cal_dvlocal_R_sparseMatrix(
        const int &current_spin,
        const double &sparse_threshold,
        LCAO_Matrix *LM,
        Parallel_Orbitals *pv,
        LCAO_Orbitals &orb,UnitCell &ucell,
        Grid_Driver &gdriver);

    // 
void process_sparse_matrix(
        int current_spin,
        double sparse_threshold,
        const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& sparseMatrix,
        LCAO_Matrix* LM,
        Parallel_Orbitals* pv,
        int dim = -1);
    private:

    //----------------------------
    // key variable 
    //----------------------------  

    // used only in vlocal.
    int spin_now = -1;


    std::vector<int> trace_lo;
    
};


inline void Gint_k::process_sparse_matrix(
    int current_spin,
    double sparse_threshold,
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& sparseMatrix,
    LCAO_Matrix* LM,
    Parallel_Orbitals* pv,
    int dim)
{
    size_t total_R_num = LM->all_R_coor.size();
    int* nonzero_num = new int[total_R_num];
    int* minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);

    size_t count = 0;
    for (const auto& R_coor : LM->all_R_coor)
    {
        auto iter = sparseMatrix.find(R_coor);
        if (iter != sparseMatrix.end())
        {
            for (const auto& row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;
        iter = sparseMatrix.find(minus_R_coor);
        if (iter != sparseMatrix.end())
        {
            for (const auto& row_loop : iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_all(minus_nonzero_num, total_R_num);

    double* tmp = new double[GlobalV::NLOCAL];
    count = 0;

    for (const auto& R_coor : LM->all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for (size_t row = 0; row < GlobalV::NLOCAL; ++row)
            {
                ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);

                auto iter = sparseMatrix.find(R_coor);
                if (iter != sparseMatrix.end())
                {
                    if (trace_lo[row] >= 0)  // 使用类成员变量 trace_lo
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (const auto& value : row_iter->second)
                            {
                                tmp[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != sparseMatrix.end())
                {
                    for (size_t col = 0; col < row; ++col)
                    {
                        if (trace_lo[col] >= 0)  // 使用类成员变量 trace_lo
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
                    for (size_t col = 0; col < GlobalV::NLOCAL; ++col)
                    {
                        if (pv->global2local_col(col) >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_threshold)
                            {
                                if (dim == -1)
                                {
                                    double& matrix_value = LM->HR_sparse[current_spin][R_coor][row][col];
                                    matrix_value += tmp[col];
                                    if (std::abs(matrix_value) <= sparse_threshold)
                                    {
                                        LM->HR_sparse[current_spin][R_coor][row].erase(col);
                                    }
                                }
                                else
                                {
                                    double* matrix_value = nullptr;
                                    if (dim == 0)
                                    {
                                        matrix_value = &LM->dHRx_sparse[current_spin][R_coor][row][col];
                                    }
                                    else if (dim == 1)
                                    {
                                        matrix_value = &LM->dHRy_sparse[current_spin][R_coor][row][col];
                                    }
                                    else if (dim == 2)
                                    {
                                        matrix_value = &LM->dHRz_sparse[current_spin][R_coor][row][col];
                                    }

                                    *matrix_value += tmp[col];
                                    if (std::abs(*matrix_value) <= sparse_threshold)
                                    {
                                        if (dim == 0)
                                        {
                                            LM->dHRx_sparse[current_spin][R_coor][row].erase(col);
                                        }
                                        else if (dim == 1)
                                        {
                                            LM->dHRy_sparse[current_spin][R_coor][row].erase(col);
                                        }
                                        else if (dim == 2)
                                        {
                                            LM->dHRz_sparse[current_spin][R_coor][row].erase(col);
                                        }
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
}


#endif
