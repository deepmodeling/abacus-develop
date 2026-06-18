#ifndef TDDEEPKS_H
#define TDDEEPKS_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_io/module_hs/cal_r_overlap_R.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#endif

#include <unordered_map>
#include <vector>

namespace hamilt
{

#ifndef __TD_DEEPKSTEMPLATE
#define __TD_DEEPKSTEMPLATE

/// The TDDeePKS class template inherits from class T.
/// It supplies the missing commutator contribution -i[V_delta, r] of the DeePKS
/// non-local potential to the rt-TDDFT velocity-gauge current (out_current=1,
/// td_stype=1). V_delta is a cross-atom projected operator
/// (V_delta_{mu,nu} = sum_iat0 <phi_mu|alpha> gedm <alpha|phi_nu>), so [V_delta,r]
/// is non-zero and contributes to the current operator. This term is otherwise
/// absent from current_term (only filled by TDEkinetic and TDNonlocal), which
/// produces a spurious DC offset in the DeePKS current/spectrum.
template <class T>
class TDDeePKS : public T
{
};

#endif

/// TDDeePKS class template specialization for OperatorLCAO<TK> base class.
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class TDDeePKS<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    TDDeePKS<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                   const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                   hamilt::HContainer<TR>* hR_in,
                                   const UnitCell* ucell_in,
                                   const LCAO_Orbitals& orb,
                                   const Grid_Driver* GridD_in
#ifdef __MLALGO
                                   ,
                                   LCAO_Deepks<TK>* ld_in
#endif
    );
    ~TDDeePKS<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() adds -i[V_delta, r] to TD_info::td_vel_op->current_term.
     * This operator does NOT modify hR (V_delta itself is added by the DeePKS
     * operator); it only forwards the shared hR_tmp pointer to the next TD
     * sub-operator (TDNonlocal) and writes the current term.
     */
    virtual void contributeHR() override;
    virtual void contributeHk(int ik) override;

    virtual void set_HR_fixed(void*) override;

  private:
    const UnitCell* ucell = nullptr;
    const LCAO_Orbitals& orb_;
    const Grid_Driver* Grid = nullptr;

#ifdef __MLALGO
    /// @brief DeePKS data handle, provides gedm and inl_index (phialpha layout)
    LCAO_Deepks<TK>* ld = nullptr;
#endif

    /// @brief shared TD real-space Hamiltonian buffer. Created by TDEkinetic and
    /// forwarded down the sub-chain. Not written here, only forwarded to the next
    /// sub-operator so that TDNonlocal accumulates Vnl into the same buffer.
    HContainer<std::complex<double>>* hR_tmp = nullptr;

    /// @brief whether the current term has been computed at least once (per instance)
    bool current_done = false;

    /// @brief search the nearest neighbor atoms (same cutoff as DeePKS V_delta:
    /// Phi[T1].Rcut + Alpha[0].Rcut, centered on the descriptor atom iat0)
    void initialize_HR(const Grid_Driver* GridD_in);

    /// @brief accumulate -i[V_delta, r] into current_term for all <I,J,R> pairs
    void calculate_current();

    /// @brief calculate the current contribution of one <I,J,R> atom pair
    void cal_current_IJR(const int& iat1,
                         const int& iat2,
                         const Parallel_Orbitals* paraV,
                         const std::vector<std::unordered_map<int, std::vector<double>>>& nlm1_all,
                         const std::vector<std::unordered_map<int, std::vector<double>>>& nlm2_all,
                         const std::vector<int>& trace_alpha_row,
                         const std::vector<int>& trace_alpha_col,
                         const std::vector<double>& gedms,
                         std::complex<double>** current_mat_p);

    /// @brief nearest neighbor atoms of each descriptor atom iat0
    std::vector<AdjacentAtomInfo> adjs_all;

    /// @brief r-weighted projection calculator <phi|alpha> and <phi|r|alpha>;
    /// the tables only depend on orbital shapes (not positions), so a single
    /// static instance, initialized once, is reused across instances and MD steps.
    static cal_r_overlap_R r_calculator;
    static bool init_done;
};

} // namespace hamilt
#endif
