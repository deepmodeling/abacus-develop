#ifndef EKINETIC_H
#define EKINETIC_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include <vector>

namespace hamilt
{

#ifndef __EKINETICTEMPLATE
#define __EKINETICTEMPLATE

/// The EKinetic class template inherits from class T
/// it is used to calculate the electronic kinetic
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T>
class EKinetic : public T
{
};

#endif

/// EKinetic class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the electronic kinetic matrix in real space and fold it to k-space
/// HR = <psi_{mu, 0}|-\Nabla^2|psi_{nu, R}>
/// HK = <psi_{mu, k}|-\Nabla^2|psi_{nu, k}> = \sum_{R} e^{ikR} HR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class EKinetic<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    /**
     * @brief Construct a new EKinetic object
     */
    EKinetic<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      HContainer<TR>* hR_in,
                                      const UnitCell* ucell_in,
                                      const std::vector<double>& orb_cutoff,
                                      const Grid_Driver* GridD_in,
                                      const TwoCenterIntegrator* intor);

    /**
     * @brief Destroy the EKinetic object
     */
    ~EKinetic<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|-\Nabla^2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    virtual void set_HR_fixed(void*) override;

    /**
     * @brief calculate force and stress for kinetic operator
     * @param cal_force whether to calculate force
     * @param cal_stress whether to calculate stress
     * @param dmR density matrix in real space
     * @param force output force matrix (nat x 3)
     * @param stress output stress matrix (3 x 3)
     */
    void cal_force_stress(const bool cal_force,
                          const bool cal_stress,
                          const HContainer<double>* dmR,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& stress);

  private:
    const UnitCell* ucell = nullptr;
    std::vector<double> orb_cutoff_;

    hamilt::HContainer<TR>* HR_fixed = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    const Grid_Driver* gridD = nullptr;

    bool allocated = false;

    bool HR_fixed_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const Grid_Driver* GridD_in);

    /**
     * @brief calculate the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * use the adjs_all to calculate the HR matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const ModuleBase::Vector3<double>& dtau,
                    TR* data_pointer);

    /**
     * @brief calculate force contribution for atom pair <I,J,R>
     */
    void cal_force_IJR(const int& iat1,
                       const int& iat2,
                       const Parallel_Orbitals* paraV,
                       const std::unordered_map<int, std::vector<double>>& nlm1_all,
                       const std::unordered_map<int, std::vector<double>>& nlm2_all,
                       const hamilt::BaseMatrix<TR>* dmR_pointer,
                       double* force1,
                       double* force2);

    /**
     * @brief calculate stress contribution for atom pair <I,J,R>
     */
    void cal_stress_IJR(const int& iat1,
                        const int& iat2,
                        const Parallel_Orbitals* paraV,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const hamilt::BaseMatrix<TR>* dmR_pointer,
                        const ModuleBase::Vector3<double>& dis1,
                        const ModuleBase::Vector3<double>& dis2,
                        double* stress);

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_all;
};

} // namespace hamilt
#endif
