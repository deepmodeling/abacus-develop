#ifndef DFTUNEW_H
#define DFTUNEW_H
#include <unordered_map>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"

namespace hamilt
{

#ifndef __DFTUNEWTEMPLATE
#define __DFTUNEWTEMPLATE

/// The DFTUNew class template inherits from class T
/// it is used to calculate the non-local pseudopotential of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T>
class DFTUNew : public T
{
};

#endif

/// DFTUNew class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the non-local pseudopotential matrix in real space and fold it to k-space
/// HR = <psi_{mu, 0}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, R}>
/// HK = <psi_{mu, k}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, k}> = \sum_{R} e^{ikR} HR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class DFTUNew<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DFTUNew<OperatorLCAO<TK, TR>>(LCAO_Matrix* LM_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* hR_in,
                                      std::vector<TK>* hK_in,
                                      const UnitCell* ucell_in,
                                      Grid_Driver* GridD_in,
                                      ModuleDFTU::DFTU* dftu_in,
                                      const Parallel_Orbitals* paraV);
    ~DFTUNew<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|beta_p1>D_{p1, p2}<beta_p2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    /// pointer of density matrix, it is a temporary implementation
    static const elecstate::DensityMatrix<TK, double>* dm_in_dftu;

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

    hamilt::HContainer<TR>* HR_fixed = nullptr;

    bool allocated = false;

    TK* HK_pointer = nullptr;

    bool HR_fixed_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    void cal_occupations(const int& iat1,
                        const int& iat2,
                        const int& T0,
                        const Parallel_Orbitals* paraV,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const double* data_pointer,
                        std::vector<double>& occupations);
    /// VU_{m, m'} = sum_{m,m'} (1/2*delta_{m, m'} - occ_{m, m'}) * U
    /// EU = sum_{m,m'} 1/2 * U * occ_{m, m'} * occ_{m', m}
    void cal_v_of_u(
      const double* occ,
      const int m_size,
      const double u_value,
      double* VU,
      double& E);

    /**
     * @brief calculate the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in HR and calculate the non-local pseudopotential matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const int& T0,
                    const Parallel_Orbitals* paraV,
                    const std::unordered_map<int, std::vector<double>>& nlm1_all,
                    const std::unordered_map<int, std::vector<double>>& nlm2_all,
                    const std::vector<double>& vu_in,
                    TR* data_pointer);

    /**
     * @brief calculate the atomic Force of <I,J,R> atom pair
     */
    void cal_force_IJR(const int& iat1,
                       const int& iat2,
                       const int& T0,
                       const Parallel_Orbitals* paraV,
                       const std::unordered_map<int, std::vector<double>>& nlm1_all,
                       const std::unordered_map<int, std::vector<double>>& nlm2_all,
                       const std::vector<double>& vu_in,
                       const hamilt::BaseMatrix<double>** dmR_pointer,
                       const int nspin,
                       double* force);
    /**
     * @brief calculate the Stress of <I,J,R> atom pair
     */
    void cal_stress_IJR(const int& iat1,
                        const int& iat2,
                        const int& T0,
                        const Parallel_Orbitals* paraV,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const std::vector<double>& vu_in,
                        const hamilt::BaseMatrix<double>** dmR_pointer,
                        const int nspin,
                        const ModuleBase::Vector3<int>& R_index1,
                        double* stress);

    std::vector<AdjacentAtomInfo> adjs_all;
};

} // namespace hamilt
#endif