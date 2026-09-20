#ifndef DFTU_NAO_OP_H
#define DFTU_NAO_OP_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h" // AdjacentAtomInfo (value member)
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_lcao/module_dftu/dftu_nao_adj.h"
#include "source_cell/module_symmetry/symm_rotation_k.h" // ModuleSymmetry::Symmetry_rotation_k (value member); LibRI-free

#include <vector>

class Plus_U_Base;
class TwoCenterIntegrator;
class UnitCell;
class K_Vectors;

namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;
} // namespace elecstate

namespace hamilt
{

/// The DFTU_onsite class template inherits from class T
/// it is used to calculate the DFT+U on-site correction of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK, TR> or OperatorPW<TK>
template <class T>
class DFTU_onsite : public T
{
};

/// DFTU_onsite class template specialization for OperatorLCAO<TK, TR> base class.
/// Adds the DFT+U on-site correction to the real-space Hamiltonian, which is
/// then folded to k-space by the OperatorLCAO machinery:
///   HR(mu,nu;I,J,R) = <phi_{mu,I,0}|chi_m> pot_onsite(m,m') <chi_m'|phi_{nu,J,R}>
///   HK = sum_R e^{ikR} HR
/// where chi_m are the Hubbard projectors of the correlated shell.
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class DFTU_onsite<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DFTU_onsite(HS_Matrix_K<TK>* hsk_in,
                const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                hamilt::HContainer<TR>* hR_in,
                const UnitCell& ucell_in,
                const Grid_Driver* gridD_in,
                const TwoCenterIntegrator* intor,
                const std::vector<double>& orb_cutoff,
                Plus_U_Base* p_dftu,
                const int nspin_in,
                const double onsite_radius,
                const elecstate::DensityMatrix<TK, double>* dm_in);
    ~DFTU_onsite() = default;

    /**
     * @brief contributeHR() calculates the HR matrix
     * <phi_{\mu, 0}|chi_m> pot_onsite(m,m') <chi_m'|phi_{\nu, R}>
     */
    void contributeHR() override;

    /// @brief provide the K_Vectors needed to restore the full-BZ density matrix
    /// under crystal symmetry (kv.kstars); set once after construction from
    /// HamiltLCAO. When left null, contributeHR() falls back to the
    /// unsymmetrized DMR (previous behavior).
    void set_kv(const K_Vectors* kv_in) { this->kv_ = kv_in; }

  private:
    const UnitCell* ucell = nullptr;

    Plus_U_Base* dftu = nullptr;

    /// @brief solver-owned density matrix providing DMR; lifetime covers each ionic step
    const elecstate::DensityMatrix<TK, double>* dm_ = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    std::vector<double> orb_cutoff_;

    /// @brief the number of spin components, 1 for no-spin, 2 for collinear spin case and 4 for non-collinear spin case
    int nspin = 0;

    /// @brief adjacent-atom lists for all Hubbard atoms; structure snapshot
    /// computed once in the constructor (operator is rebuilt every ionic step)
    std::vector<AdjacentAtomInfo> adjs_all;
    /// @brief cached <phi|alpha^I> overlap values; structure snapshot computed
    /// once in the constructor, reused across SCF iterations of one ionic step
    DFTU_LCAO::NlmTot nlm_tot;

    /// @brief K_Vectors of the current run; only used (when non-null) to restore
    /// the full-BZ density matrix under crystal symmetry before it is read into
    /// the DFT+U occupation matrix. Set via set_kv() after construction.
    const K_Vectors* kv_ = nullptr;

    /// @brief LibRI-independent k-space rotation machinery (shared base of EXX's own
    /// Symmetry_rotation) to reconstruct D(k) at every k-star member from D(k_ibz);
    /// built lazily (once) the first time it is needed (rebuilt every ionic step).
    ModuleSymmetry::Symmetry_rotation_k symrot_;
    bool symrot_built_ = false;
};

} // namespace hamilt
#endif
