#ifndef DFTPLUSU_H
#define DFTPLUSU_H
#include <unordered_map>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "dftu.hpp"

namespace hamilt
{

/// DFTU class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the non-local pseudopotential matrix in real space and fold it to k-space
/// HR = <psi_{mu, 0}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, R}>
/// HK = <psi_{mu, k}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, k}> = \sum_{R} e^{ikR} HR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class DFTU<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DFTU<OperatorLCAO<TK, TR>>(LCAO_Matrix* lm_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* hR_in,
                                      std::vector<TK>* hK_in,
                                      const UnitCell& ucell_in,
                                      Grid_Driver* gridD_in,
                                      ModuleDFTU::DFTU* dftu_in,
                                      const Parallel_Orbitals& paraV);
    ~DFTU<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|beta_p1>D_{p1, p2}<beta_p2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    /// calculate force and stress for DFT+U 
    void cal_force_stress(
      const bool cal_force, 
      const bool cal_stress, 
      ModuleBase::matrix& force, 
      ModuleBase::matrix& stress);

  private:
    const UnitCell* ucell = nullptr;

    ModuleDFTU::DFTU* dftu = nullptr;

    hamilt::HContainer<TR>* HR = nullptr;

    TK* HK_pointer = nullptr;

    /// @brief the number of spin components, 1 for no-spin, 2 for collinear spin case and 4 for non-collinear spin case
    int nspin = 0;
    /// @brief the current spin index for nspin==2 to calculate spin-up and spin-down separately
    int current_spin = 0;

    /**
     * @brief search the nearest neighbor atoms and save them into this->adjs_all
     * the size of HR will not change in DFTU, 
     * because I don't want to expand HR larger than Nonlocal operator caused by DFTU
     */
    void initialize_HR(Grid_Driver* gridD_in, const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the <phi|alpha^I> overlap values and save them in this->nlm_tot
     * it will be reused in the calculation of calculate_HR()
    */
    void cal_nlm_all(const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the occ_mm' = \sum_R DMR*<phi_0|alpha^I_m'><alpha^I_m'|phi_R> matrix for each atom to add U
    */
    void cal_occ(const int& iat1,
                const int& iat2,
                const Parallel_Orbitals* paraV,
                const std::unordered_map<int, std::vector<double>>& nlm1_all,
                const std::unordered_map<int, std::vector<double>>& nlm2_all,
                const double* data_pointer,
                std::vector<double>& occupations);

    /// transfer VU format from pauli matrix to normal for non-collinear spin case
    void transfer_vu(std::vector<double>& vu_tmp, std::vector<TR>& vu);
    /// VU_{m, m'} = sum_{m,m'} (1/2*delta_{m, m'} - occ_{m, m'}) * U
    /// EU = sum_{m,m'} 1/2 * U * occ_{m, m'} * occ_{m', m}
    void cal_v_of_u(
      const std::vector<double>& occ,
      const int m_size,
      const double u_value,
      double* vu,
      double& eu);

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const std::unordered_map<int, std::vector<double>>& nlm1_all,
                    const std::unordered_map<int, std::vector<double>>& nlm2_all,
                    const std::vector<TR>& vu_in,
                    TR* data_pointer);

    /**
     * @brief calculate the atomic Force of <I,J,R> atom pair
     */
    void cal_force_IJR(const int& iat1,
                       const int& iat2,
                       const Parallel_Orbitals* paraV,
                       const std::unordered_map<int, std::vector<double>>& nlm1_all,
                       const std::unordered_map<int, std::vector<double>>& nlm2_all,
                       const std::vector<double>& vu_in,
                       const hamilt::BaseMatrix<double>** dmR_pointer,
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
                        const std::vector<double>& vu_in,
                        const hamilt::BaseMatrix<double>** dmR_pointer,
                        const int nspin,
                        const ModuleBase::Vector3<double>& dis1,
                        const ModuleBase::Vector3<double>& dis2,
                        double* stress);

    std::vector<AdjacentAtomInfo> adjs_all;
    /// @brief if the nlm_tot is calculated
    bool precal_nlm_done = false;
    /// @brief the overlap values for all [atoms][nerghbors][orb_index(iw) in NAOs][m of target_l in Projectors]
    std::vector<std::vector<std::unordered_map<int, std::vector<double>>>> nlm_tot;
};

} // namespace hamilt
#endif