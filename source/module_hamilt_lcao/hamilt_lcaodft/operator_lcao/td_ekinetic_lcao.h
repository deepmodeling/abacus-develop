#ifndef TDEKINETIC_H
#define TDEKINETIC_H
#include "module_base/timer.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/klist.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#include "operator_lcao.h"
#include <vector>


namespace hamilt
{

#ifndef __TDEKINETICTEMPLATE
#define __TDEKINETICTEMPLATE

/// The TDEkinetic class template inherits from class T
/// It is used to calculate correction term of kinetic energy in time-dependent DFT
template <class T>
class TDEkinetic : public T
{
};

#endif

/// TDEkinetic class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate correction term of kinetic energy in time-dependent DFT
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian

template <typename TK, typename TR>
class TDEkinetic<OperatorLCAO<TK,TR>> : public OperatorLCAO<TK, TR>
{
  public:
    TDEkinetic<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                     hamilt::HContainer<TR>* hR_in,
                                     const K_Vectors* kv_in,
                                     const UnitCell* ucell_in,
                                     const std::vector<double>& orb_cutoff,
                                     const Grid_Driver* GridD_in,
                                     const TwoCenterIntegrator* intor);
    ~TDEkinetic();

    virtual void contributeHR() override;
    virtual void contributeHk(int ik) override;
    /// @brief init two center integrals and vector potential for td_ekintic term
    void init_td();

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const Grid_Driver* GridD);

    /**
     * @brief initialize HR_tmp
     * Allocate the memory for HR_tmp with the same size as HR
    */
    void initialize_HR_tmp();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const ModuleBase::Vector3<double>& dtau,
                    std::complex<double>* hr_mat_p,
                    std::complex<double>** current_mat_p);

    /**
     * @brief calculate the ekinetic matrix correction term in tddft with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in HR and calculate the ekinetic matrix
     */
    void calculate_HR();
    virtual void set_HR_fixed(void*) override;

    TD_Velocity td_velocity;

  private:
    const UnitCell* ucell = nullptr;
    std::vector<double> orb_cutoff_;
    
    HContainer<TR>* SR = nullptr;
    /// @brief Store real space hamiltonian. TD term should include imaginary part, thus it has to be complex type. Only shared between TD operators.
    HContainer<std::complex<double>>* hR_tmp = nullptr;
    const Grid_Driver* Grid = nullptr;

    const K_Vectors* kv;
    /// @brief correction term iA⋅∇
    void td_ekinetic_scalar(std::complex<double>* Hloc, const TR& Sloc, int nnr);
    /// @brief correction term A^2*S
    void td_ekinetic_grad(std::complex<double>* Hloc, int nnr, ModuleBase::Vector3<double> grad_overlap);

    const TwoCenterIntegrator* intor_ = nullptr;

    /// @brief Store the vector potential for td_ekinetic term
    ModuleBase::Vector3<double> cart_At;

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_all;

    bool hR_tmp_done = false;
    bool allocated = false;
    bool output_hR_done = false;
};

} // namespace hamilt
#endif
