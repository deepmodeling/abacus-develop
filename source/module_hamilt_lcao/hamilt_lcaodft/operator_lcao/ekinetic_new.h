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

/// The OverlapNew class template inherits from class T
/// it is used to calculate the overlap of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T, typename TR>
class EkineticNew : public T
{
};

#endif

/// OverlapNew class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the overlap matrix in real space and fold it to k-space
/// HR = <psi_{mu, 0}|-\Nabla^2|psi_{nu, R}>
/// HK = <psi_{mu, k}|-\Nabla^2|psi_{nu, k}> = \sum_{R} e^{ikR} HR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class EkineticNew<OperatorLCAO<TK>, TR> : public OperatorLCAO<TK>
{
  public:
    EkineticNew<OperatorLCAO<TK>, TR>(LCAO_Matrix* LM_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* SR_in,
                                      TK* HK_pointer_in,
                                      const UnitCell* ucell_in,
                                      Grid_Driver* GridD_in,
                                      const Parallel_Orbitals* paraV);

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    const UnitCell* ucell = nullptr;

    hamilt::HContainer<TR>* HR = nullptr;

    TK* HK_pointer = nullptr;

    bool HR_fixed_done = false;

    /**
     * @brief initialize SR, search the nearest neighbor atoms
     * HContainer is used to store the overlap matrix with specific <I,J,R> atom-pairs
     * the size of SR will be fixed after initialization
     */
    void initialize_HR(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the overlap matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in SR and calculate the overlap matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the SR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const ModuleBase::Vector3<double>& dtau,
                    TR* data_pointer);
};

} // namespace hamilt
#endif