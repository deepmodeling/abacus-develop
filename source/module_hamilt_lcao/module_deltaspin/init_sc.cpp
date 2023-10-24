#include "spin_constrain.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::set_solver_parameters(
    int nspin_in,
    K_Vectors kv_in,
    hsolver::HSolver<std::complex<double>, psi::DEVICE_CPU>* phsol_in,
    hamilt::Hamilt<std::complex<double>, psi::DEVICE_CPU>* p_hamilt_in,
    psi::Psi<std::complex<double>>* psi_in,
    elecstate::ElecState* pelec_in,
    std::string KS_SOLVER_in,
    LCAO_Matrix* LM_in)
{
    /// set nspin
    this->set_nspin(nspin_in);
    this->kv_ = kv_in;
    this->phsol = phsol_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
    this->KS_SOLVER = KS_SOLVER_in;
    this->LM = LM_in;
}

// init sc
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::init_sc(const UnitCell& ucell,
                                            int NPOL,
                                            std::string sc_file,
                                            Parallel_Orbitals* ParaV_in,
                                            int nspin_in,
                                            double sc_thr_in,
                                            int nsc_in,
                                            int nsc_min_in,
                                            double alpha_trial_in,
                                            double sccut_in,
                                            bool decay_grad_switch_in,
                                            K_Vectors kv_in,
                                            std::string KS_SOLVER_in,
                                            LCAO_Matrix* LM_in,
                                            hsolver::HSolver<std::complex<double>, psi::DEVICE_CPU>* phsol_in,
                                            hamilt::Hamilt<std::complex<double>, psi::DEVICE_CPU>* p_hamilt_in,
                                            psi::Psi<std::complex<double>>* psi_in,
                                            elecstate::ElecState* pelec_in)
{
    this->set_input_parameters(sc_thr_in, nsc_in, nsc_min_in, alpha_trial_in, sccut_in, decay_grad_switch_in);
    this->set_ParaV(ParaV_in);
    this->set_solver_parameters(nspin_in, kv_in, phsol_in, p_hamilt_in, psi_in, pelec_in, KS_SOLVER_in, LM_in);
    this->bcast_ScData(ucell, NPOL, sc_file);
}