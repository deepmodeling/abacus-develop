#ifndef DELTA_SPIN_LCAO_H
#define DELTA_SPIN_LCAO_H

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include <unordered_map>

namespace hamilt
{

#ifndef __DELTASPINTEMPLATE
#define __DELTASPINTEMPLATE

template <class T>
class DeltaSpin : public T
{
};

#endif

template <typename TK, typename TR>
class DeltaSpin<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DeltaSpin<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* hR_in,
                                      const UnitCell& ucell_in,
                                      Grid_Driver* gridD_in,
                                      const TwoCenterIntegrator* intor,
                                      const std::vector<double>& orb_cutoff);
    ~DeltaSpin<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|beta_p1>D_{p1, p2}<beta_p2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    /**
     * @brief calculate the magnetization moment for each atom
     * @param dmR the density matrix in real space
     * @return the magnetization moment for each atom
    */
    std::vector<double> cal_moment(const HContainer<double>* dmR, const std::vector<ModuleBase::Vector3<int>>& constrain);

    /**
     * @brief set the update_lambda_ to true, which means the lambda will be updated in the next contributeHR()
    */
    void update_lambda()
    {
        for(int is=0;is<this->spin_num;is++)
        {
            this->update_lambda_[is] = true;
        }
    }

    /// calculate force and stress for DFT+U
    void cal_force_stress(const bool cal_force,
                          const bool cal_stress,
                          const HContainer<double>* dmR,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& stress);

  private:
    const UnitCell* ucell = nullptr;

    Grid_Driver* gridD = nullptr;

    const Parallel_Orbitals* paraV = nullptr;

    hamilt::HContainer<TR>* HR = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    std::vector<double> orb_cutoff_;

    /// @brief the number of spin components, 1 for no-spin, 2 for collinear spin case and 4 for non-collinear spin case
    int nspin = 0;

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const std::unordered_map<int, std::vector<double>>& nlm1_all,
                    const std::unordered_map<int, std::vector<double>>& nlm2_all,
                    TR* data_pointer);

    /**
     * @brief calculate the prepare HR for each atom
     * pre_hr^I = \sum_{lm}<phi_mu|alpha^I_{lm}><alpha^I_{lm}|phi_{nu,R}>
     */
    void cal_pre_HR();

    /**
     * @brief calculate the constaint atom list
    */
    void cal_constraint_atom_list(const std::vector<ModuleBase::Vector3<int>>& constraints);

    /**
     * @brief calculate the atomic magnetization moment for each <IJR>
    */
    void cal_moment_IJR(const double* dmR, 
                        const TR* hr, 
                        const int row_size,
                        const int col_size,
                        double* moment);

    /**
     * @brief calculate the atomic Force of <I,J,R> atom pair
     */
    void cal_force_IJR(const int& iat1,
                       const int& iat2,
                       const Parallel_Orbitals* paraV,
                       const std::unordered_map<int, std::vector<double>>& nlm1_all,
                       const std::unordered_map<int, std::vector<double>>& nlm2_all,
                       const hamilt::BaseMatrix<double>* dmR_pointer,
                       const ModuleBase::Vector3<double>& lambda,
                       const int nspin,
                       double* force1,
                       double* force2);
    /**
     * @brief calculate the Stress of <I,J,R> atom pair
     */
    void cal_stress_IJR(const int& iat1,
                        const int& iat2,
                        const Parallel_Orbitals* paraV,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const hamilt::BaseMatrix<double>* dmR_pointer,
                        const ModuleBase::Vector3<double>& lambda,
                        const int nspin,
                        const ModuleBase::Vector3<double>& dis1,
                        const ModuleBase::Vector3<double>& dis2,
                        double* stress);

    /**
     * @brief calculate the array of coefficient of lambda * d\rho^p/drho^{\sigma\sigma'}
    */
    void pre_coeff_array(const std::vector<TR>& coeff, const int row_size, const int col_size);

    std::vector<bool> constraint_atom_list;
    std::vector<hamilt::HContainer<TR>*> pre_hr;

    std::vector<double> tmp_dmr_memory;
    std::vector<TR> tmp_coeff_array;
    std::vector<double> lambda_save;

    bool initialized = false;
    int spin_num = 1;
    std::vector<bool> update_lambda_;
};

}

#endif