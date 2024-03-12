#ifndef OVERLAPNEW_H
#define OVERLAPNEW_H
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{

#ifndef __OVERLAPNEWTEMPLATE
#define __OVERLAPNEWTEMPLATE

/// The OverlapNew class template inherits from class T
/// it is used to calculate the overlap of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T>
class OverlapNew : public T
{
};

#endif

/// OverlapNew class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the overlap matrix in real space and fold it to k-space
/// SR = <psi_{mu, 0}|psi_{nu, R}>
/// SK = <psi_{mu, k}|psi_{nu, k}> = \sum_{R} e^{ikR} SR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class OverlapNew<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    OverlapNew<OperatorLCAO<TK, TR>>(LCAO_Matrix* LM_in,
                                     const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                     hamilt::HContainer<TR>* hR_in,
                                     std::vector<TK>* hK_in,
                                     hamilt::HContainer<TR>* SR_in,
                                     std::vector<TK>* SK_pointer_in,
                                     const UnitCell* ucell_in,
                                     Grid_Driver* GridD_in,
                                     const Parallel_Orbitals* paraV);

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

    TK* getSk();

  private:
    const UnitCell* ucell = nullptr;

    hamilt::HContainer<TR>* SR = nullptr;

    std::vector<TK>* SK_pointer = nullptr;

    bool SR_fixed_done = false;

    /**
     * @brief initialize SR, search the nearest neighbor atoms
     * HContainer is used to store the overlap matrix with specific <I,J,R> atom-pairs
     * the size of SR will be fixed after initialization
     */
    void initialize_SR(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the overlap matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in SR and calculate the overlap matrix
     */
    void calculate_SR();

    /**
     * @brief calculate the SR local matrix of <I,J,R> atom pair
     */
    void cal_SR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const ModuleBase::Vector3<double>& dtau,
                    TR* data_pointer);

    // if k vector is not changed, then do nothing and return
    // default of kvec_d_old is (-10,-10,-10), which is not a valid k vector
    ModuleBase::Vector3<double> kvec_d_old = ModuleBase::Vector3<double>(-10,-10,-10);
};

} // namespace hamilt
#endif