#ifndef DELTA_SPIN_LCAO_H
#define DELTA_SPIN_LCAO_H

#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include <unordered_map>
#include <complex>

class Parallel_Blkcyclic;

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
                                    const Grid_Driver* gridD_in,
                                    const TwoCenterIntegrator* intor,
                                    const std::vector<double>& orb_cutoff);

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

    /// @brief Reset initialization state to allow re-constraint with new constrain array
    void reset_initialized()
    {
        this->initialized = false;
    }

    /**
     * @brief set the update_lambda_ to true, which means the lambda will be updated in the next contributeHR()
    */
    void update_lambda()
    {
        for(int is=0;is<this->spin_num;is++)
        {
            this->update_lambda_[is] = true;
        }
        // Reset sc_hr_done so contributeHR() recalculates DeltaSpin HR
        // in the next k-point loop (avoids accumulation across k-points)
        this->sc_hr_done = false;
    }

    /**
     * @brief Shadow set_current_spin to reset sc_hr_done on spin switch (nspin=2).
     * In the lambda loop, refresh_times=0 so the shared hr_done is NOT reset on
     * spin switch. sc_hr_done must be reset here so each spin's HR is computed
     * independently.
     */
    void set_current_spin(const int current_spin_in) override
    {
        if (this->current_spin != current_spin_in)
        {
            this->sc_hr_done = false;
        }
        OperatorLCAO<TK, TR>::set_current_spin(current_spin_in);
    }

    /// calculate force and stress for DFT+U
    void cal_force_stress(const bool cal_force,
                          const bool cal_stress,
                          const HContainer<double>* dmR,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& stress);

    /// Compute P_I_sub(k) = D_I(k)^dag D_I(k) for all constrained atoms.
    void cal_PI_sub(const Parallel_Blkcyclic& matrix_layout,
                    const Parallel_Blkcyclic& wfc_layout,
                    const ModuleBase::Vector3<double>& kvec_d,
                    const std::complex<double>* psi_k,
                    const int nbands_global,
                    std::vector<std::vector<std::complex<double>>>& PI_sub) const;

  private:
    const UnitCell* ucell = nullptr;

    const Grid_Driver* gridD = nullptr;

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
                       const HContainer<double>& dmR,
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
                        const HContainer<double>& dmR,
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
    std::vector<hamilt::HContainer<TR>> pre_hr;

    std::vector<double> tmp_dmr_memory;
    std::vector<TR> tmp_coeff_array;
    std::vector<double> lambda_save;

    bool initialized = false;
    int spin_num = 1;
    std::vector<bool> update_lambda_;
    /// Independent HR completion flag for DeltaSpin, decoupled from
    /// the shared OperatorLCAO::hr_done to avoid cross-k-point accumulation.
    bool sc_hr_done = false;

    /// Saved alpha-orbital overlaps used by cal_PI_sub.
    struct BI_AdjacentData {
        int iat_adj;
        ModuleBase::Vector3<int> R_index;
        std::unordered_map<int, std::vector<double>> nlm;
    };
    std::vector<std::vector<BI_AdjacentData>> B_I_data;
    std::vector<int> B_I_nproj;
};

}

#endif
