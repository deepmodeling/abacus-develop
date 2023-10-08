#ifndef GINT_K_H
#define GINT_K_H

#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "grid_technique.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_elecstate/module_charge/charge.h"
#include "gint.h"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "module_base/abfs-vector3_order.h"

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
    void folding_vl_k(const int &ik, LCAO_Matrix* LM, const std::vector<ModuleBase::Vector3<double>>& kvec_d);

    /**
     * @brief transfer pvpR to this->hRGint
     * then pass this->hRGint to Veff<OperatorLCAO>::hR
    */
    void transfer_pvpR(hamilt::HContainer<double> *hR);
    void transfer_pvpR(hamilt::HContainer<std::complex<double>> *hR);

    //------------------------------------------------------
    // in gint_k_env.cpp 
    //------------------------------------------------------
    // calculate the envelop function via grid integrals
    void cal_env_k(int ik,
                   const std::complex<double>* psi_k,
                   double* rho,
                   const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                   const std::vector<ModuleBase::Vector3<double>>& kvec_d);

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
        LCAO_Matrix *LM);

    void distribute_pvpR_soc_sparseMatrix(
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t,
        std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix,
        LCAO_Matrix *LM);

    void cal_vlocal_R_sparseMatrix(
        const int &current_spin,
        const double &sparse_threshold,
        LCAO_Matrix *LM);

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
        LCAO_Matrix *LM);

    void distribute_pvdpR_soc_sparseMatrix(
        const int dim,
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t, std::map<size_t, std::complex<double>>>> &pvdpR_soc_sparseMatrix,
        LCAO_Matrix *LM);

    void cal_dvlocal_R_sparseMatrix(
        const int &current_spin,
        const double &sparse_threshold,
        LCAO_Matrix *LM);

    private:

    //----------------------------
    // key variable 
    //----------------------------  

    // used only in vlocal.
    int spin_now = -1;
    
};

#endif
