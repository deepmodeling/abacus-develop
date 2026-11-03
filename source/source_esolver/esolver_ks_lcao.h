#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H

#include "esolver_ks.h"
#include "source_basis/module_nao/two_center_bundle.h" // nao basis
#include "source_cell/record_adj.h"                    // adjacent atoms
#include "source_hamilt/module_gint/gint_info.h"
#include "source_hamilt/module_gint/gint_prec_ctrl.h"
#include "source_hamilt/module_xc/exx_info.h" // LCAO owns full Exx_Info
#include "source_lcao/allocate_dm.h"          // mohan add 2025-10-30
#include "source_lcao/module_rdmft/rdmft.h"   // rdmft
#include "source_lcao/module_rt/td_info.h"
#include "source_lcao/setup_deepks.h" // for deepks, mohan add 20251008
#include "source_lcao/setup_exx.h"    // for exx, mohan add 20251008

#include <complex>
#include <memory>

//-----------------------------------
// ESolver for LCAO
//-----------------------------------
namespace ModuleESolver
{
//! Forward declaration for the linear-response solver.
template <typename T, typename TR>
class ESolver_LR;

//! Kohn-Sham energy solver using the LCAO basis.
template <typename TK, typename TR>
class ESolver_KS_LCAO : public ESolver_KS
{
  public:
    //! Initialize the solver identity and DFT+U storage.
    ESolver_KS_LCAO();
    //! Release the electronic wave functions.
    ~ESolver_KS_LCAO();

    //! Initialize the LCAO basis, electronic state, and potentials before all steps.
    void before_all_runners(BaseCell& basecell, const Input_para& inp) override;

    //! Return the electronic total energy in Ry.
    double cal_energy() override;

    //! Calculate the requested forces and stress, caching the stress for cal_stress().
    void cal_force(BaseCell& basecell, ModuleBase::matrix& force) override;

    //! Return the cached stress, calling cal_force() first when needed.
    void cal_stress(BaseCell& basecell, ModuleBase::matrix& stress) override;

    //! Perform final LCAO output and postprocessing after all steps.
    void after_all_runners(BaseCell& basecell) override;

  protected:
    //! Optional solver-owned RT state passed to LCAO consumers.
    virtual TD_info* get_td_info() const { return nullptr; }

    //! Prepare an ionic step using istep for both initialization and initial-state output.
    virtual void before_scf(UnitCell& ucell, const int istep) override;

    //! Prepare an ionic step with a separate index for initial charge and potential output.
    void before_scf(UnitCell& ucell, const int istep, const int output_step);

    //! Prepare charge mixing and Hamiltonian updates for an SCF iteration.
    virtual void iter_init(UnitCell& ucell, const int istep, const int iter) override;

    //! Solve the Hamiltonian and update the charge density and energy terms.
    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    //! Finish an SCF iteration with energy corrections, mixing, convergence checks, and output.
    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    //! Write the SCF results and release adjacency records when forces and stress are not requested.
    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    //! Run specialized LCAO calculations outside the standard SCF workflow.
    virtual void others(BaseCell& basecell, const int istep) override;

    //! Electronic wave functions (moved from base class)
    psi::Psi<TK>* psi = nullptr;

    //! Store information about adjacent atoms for the distributed matrix layout.
    Record_adj RA;

    //! Search for adjacent atoms and provide their positions and cell indices.
    Grid_Driver gd;

    //! NAO orbitals: 2D block-cyclic distribution info
    Parallel_Orbitals pv;

    //! GintInfo: used to store basic grid-integration information for module_gint
    std::unique_ptr<ModuleGint::GintInfo> gint_info_;

    //! Store NAO orbitals and their radial tables.
    LCAO_Orbitals orb_;

    //! NAO orbitals: two-center integrations
    TwoCenterBundle two_center_bundle_;

    //! Own the density matrix, mohan add 2025-10-30.
    module_dm::Setup_DM<TK> dmat;

    //! For DeePKS calculations, mohan add 2025-10-08.
    Setup_DeePKS<TK> deepks;

    //! Full EXX info for LCAO (includes info_ri, info_opt_abfs, info_lip)
    Exx_Info exx_info_;

    //! For exact-exchange calculations, mohan add 2025-10-08.
    Exx_NAO<TK> exx_nao;

    //! For RDMFT calculations, added by jghan, 2024-03-16
    rdmft::RDMFT<TK, TR> rdmft_solver;

    //! For linear-response TDDFT
    friend class ESolver_LR<double, double>;
    //! For linear-response TDDFT
    friend class ESolver_LR<std::complex<double>, double>;

    //! Cache the stress computed by cal_force() for the separate cal_stress() interface.
    ModuleBase::matrix scs;
    //! Track whether cal_force() has run since the last cal_stress() call.
    bool have_force = false;

    //! Control grid-integration precision and its switching during SCF iterations.
    GintPrecisionController gint_precision_controller_;
};
} // namespace ModuleESolver
#endif
