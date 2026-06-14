// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "esolver_dfpt_pw.h"
#include "source_base/tool_quit.h"

namespace ModuleESolver {

ESolver_DFPT_PW::ESolver_DFPT_PW() {
    this->classname = "ESolver_DFPT_PW";
    this->basisname = "PW";
    gs_done_ = false;
    dfpt_ = nullptr;
}

ESolver_DFPT_PW::~ESolver_DFPT_PW() {
    if (dfpt_ != nullptr) {
        delete dfpt_;
        dfpt_ = nullptr;
    }
}

void ESolver_DFPT_PW::before_all_runners(UnitCell& ucell, const Input_para& inp) {
    ModuleBase::TITLE("ESolver_DFPT_PW", "before_all_runners");
    
    ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::before_all_runners(ucell, inp);
    
    init_dfpt(ucell);
}

void ESolver_DFPT_PW::runner(UnitCell& ucell, const int istep) {
    ModuleBase::TITLE("ESolver_DFPT_PW", "runner");
    
    if (!gs_done_) {
        run_gs(ucell);
        gs_done_ = true;
    }
    
    if (dfpt_ != nullptr) {
        dfpt_->run();
    }
    
    run_post_process(ucell);
}

void ESolver_DFPT_PW::after_all_runners(UnitCell& ucell) {
    ModuleBase::TITLE("ESolver_DFPT_PW", "after_all_runners");
    
    ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::after_all_runners(ucell);
}

void ESolver_DFPT_PW::run_gs(UnitCell& ucell) {
    ModuleBase::TITLE("ESolver_DFPT_PW", "run_gs");
    
    ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::runner(ucell, 0);
}

void ESolver_DFPT_PW::init_dfpt(UnitCell& ucell) {
    ModuleBase::TITLE("ESolver_DFPT_PW", "init_dfpt");
    
    dfpt_ = new ModuleDFPT::DFPT_PW();
    
//    dfpt_->init(ucell, *this->stp.psi, this->pelec->nelec, PARAM.inp.ecutwfc);
    
    dfpt_->set_parameters("dfpt.in");
    
    dfpt_->set_qmesh(1, 1, 1);
    
    dfpt_->set_conv_thr(1e-8);
    dfpt_->set_max_iter(100);
}

void ESolver_DFPT_PW::run_post_process(UnitCell& ucell) {
    ModuleBase::TITLE("ESolver_DFPT_PW", "run_post_process");
}

} // namespace ModuleESolver
