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

namespace ModuleESolver {

class ESolver_DFPT_PW : public ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU> {
public:
    ESolver_DFPT_PW();
    ~ESolver_DFPT_PW();
    
    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;
    void runner(UnitCell& ucell, const int istep) override;
    void after_all_runners(UnitCell& ucell) override;
    
protected:
    ModuleDFPT::DFPT_PW* dfpt_ = nullptr;
    
    bool gs_done_ = false;
    
    void run_gs(UnitCell& ucell);
    
    void init_dfpt(UnitCell& ucell);
    
    void run_post_process(UnitCell& ucell);
};

} // namespace ModuleESolver

#endif // ESOLVER_DFPT_PW_H
