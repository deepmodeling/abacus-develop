#ifndef TDNONLOCAL_H
#define TDNONLOCAL_H
#include <unordered_map>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"

namespace hamilt
{

#ifndef __TD_NONLOCALTEMPLATE
#define __TD_NONLOCALTEMPLATE

/// The TDNonlocal class template inherits from class T
/// It is used to calculate correction term of non-local pseudopotential in time-dependent DFT
template <class T>
class TDNonlocal : public T
{
};

#endif

/// TDNonlocal class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate correction term of non-local pseudopotential in time-dependent DFT
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class TDNonlocal<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    TDNonlocal<OperatorLCAO<TK, TR>>(LCAO_Matrix* LM_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* hR_in,
                                      std::vector<TK>* hK_in,
                                      const UnitCell* ucell_in,
                                      Grid_Driver* GridD_in,
                                      const Parallel_Orbitals* paraV);
    ~TDNonlocal<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|beta_p1>D_{p1, p2}<beta_p2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;
    virtual void contributeHk(int ik) override;

    virtual void set_HR_fixed(void*) override;

  private:
    const UnitCell* ucell = nullptr;

    HContainer<TR>* HR = nullptr;
    /// @brief Store real space hamiltonian. TD term should include imaginary part, thus it has to be complex type. Only shared between TD operators.
    HContainer<std::complex<double>>* hR_tmp = nullptr;
    Grid_Driver* Grid = nullptr;

    bool allocated = false;

    TK* HK_pointer = nullptr;

    bool hR_tmp_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    /**
     * @brief initialize HR_tmp
     * Allocate the memory for HR_tmp with the same size as HR
    */
    void initialize_HR_tmp(const Parallel_Orbitals* paraV);
    /// @brief init vector potential for td_nonlocal term
    void init_td(void);
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
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm1_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm2_all,
                    std::complex<double>* data_pointer);

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_all;
    /// @brief Store the vector potential for td_nonlocal term
    ModuleBase::Vector3<double> cart_At;
};

} // namespace hamilt
#endif
