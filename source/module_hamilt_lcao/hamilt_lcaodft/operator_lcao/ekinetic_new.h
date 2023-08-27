#ifndef EKINETICNEW_H
#define EKINETICNEW_H
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{

#ifndef __EKINETICNEWTEMPLATE
#define __EKINETICNEWTEMPLATE

/// The EkineticNew class template inherits from class T
/// it is used to calculate the electronic kinetic
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T, typename TR>
class EkineticNew : public T
{
};

#endif

/// EkineticNew class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the electronic kinetic matrix in real space and fold it to k-space
/// HR = <psi_{mu, 0}|-\Nabla^2|psi_{nu, R}>
/// HK = <psi_{mu, k}|-\Nabla^2|psi_{nu, k}> = \sum_{R} e^{ikR} HR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class EkineticNew<OperatorLCAO<TK>, TR> : public OperatorLCAO<TK>
{
  public:
    /**
     * @brief Construct a new EkineticNew object
     */
    EkineticNew<OperatorLCAO<TK>, TR>(LCAO_Matrix* LM_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* HR_in,
                                      TK* HK_pointer_in,
                                      const UnitCell* ucell_in,
                                      Grid_Driver* GridD_in,
                                      const Parallel_Orbitals* paraV);

    /**
     * @brief Destroy the EkineticNew object
     */
    ~EkineticNew<OperatorLCAO<TK>, TR>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|-\Nabla^2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    /**
     * @brief contributeHk() is used to calculate the HK matrix
     * <phi_{\mu, k}|-\Nabla^2|phi_{\nu, k}> = \sum_{R} e^{ikR} HR
     * in this class, TR has been added into total HR, skip this step.
     * @param ik: index of k-point
     */
    virtual void contributeHk(int ik) override;

    virtual void set_HR_fixed(void*) override;

  private:
    const UnitCell* ucell = nullptr;

    hamilt::HContainer<TR>* HR = nullptr;

    hamilt::HContainer<TR>* HR_fixed = nullptr;

    bool allocated = false;

    TK* HK_pointer = nullptr;

    bool HR_fixed_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in HR and calculate the electronic kinetic matrix
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

    std::vector<AdjacentAtomInfo> adjs_all;
};

} // namespace hamilt
#endif