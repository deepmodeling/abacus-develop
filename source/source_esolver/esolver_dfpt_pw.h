// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef ESOLVER_DFPT_PW_H
#define ESOLVER_DFPT_PW_H

#include "esolver_ks_pw.h"
#include "source_pw/module_dfpt/dfpt_pw.h"

namespace ModuleDFPT {
class XC_First_Order;
}

namespace ModuleESolver
{

class ESolver_DFPT_PW : public ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>
{
  public:
    ESolver_DFPT_PW();
    ~ESolver_DFPT_PW();

    void before_all_runners(BaseCell& basecell, const Input_para& inp) override;
    void runner(BaseCell& basecell, const int istep) override;
    void after_all_runners(BaseCell& basecell) override;

  protected:
    ModuleDFPT::DFPT_PW* dfpt_ = nullptr;

    ///< first-order XC kernel adapter over elecstate::PotXC_FDM (C7),
    ///< owned here so module_dfpt stays free of pot_xc_fdm.h dependencies
    ModuleDFPT::XC_First_Order* xc_adapter_ = nullptr;

    bool gs_done_ = false;

    bool dfpt_wired_ = false;

    ///< ground-state scalars captured from Input_para in before_all_runners
    ///< (rule 1: passed explicitly instead of re-reading the global record
    ///< in init_dfpt)
    int nspin_ = 1;

    double nelec_ = 0.0;

    double ecutwfc_ = 0.0;

    bool dft_plus_u_ = false;

    void run_gs(UnitCell& ucell);

    /// wires DFPT_PW with the converged ground state; called after run_gs
    /// (the injected veff/XC reference data only exist once the GS SCF is
    /// done)
    void init_dfpt(UnitCell& ucell);

    void run_post_process(UnitCell& ucell);
};

} // namespace ModuleESolver

#endif // ESOLVER_DFPT_PW_H
