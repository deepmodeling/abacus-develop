#include "spin_constrain.h"

// init sc
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::init_sc(
    double sc_thr_in,
    int nsc_in,
    int nsc_min_in,
    double alpha_trial_in,
    double sccut_in,
    bool decay_grad_switch_in,
    const UnitCell& ucell,
    std::string sc_file,
    int NPOL,
    Parallel_Orbitals* ParaV_in,
    int nspin_in,
    K_Vectors kv_in,
    std::string KS_SOLVER_in,
    LCAO_Matrix* LM_in,
    hsolver::HSolver<std::complex<double>, psi::DEVICE_CPU>* phsol_in,
    hamilt::Hamilt<std::complex<double>, psi::DEVICE_CPU>* p_hamilt_in,
    psi::Psi<std::complex<double>>* psi_in,
    elecstate::ElecState* pelec_in)
{
    this->set_input_parameters(sc_thr_in, nsc_in, nsc_min_in, alpha_trial_in, sccut_in, decay_grad_switch_in);
    std::map<int, int> atomCounts = ucell.get_atomCounts();
    std::map<int, int> orbitalCounts = ucell.get_orbitalCounts();
    this->set_orb_counts(atomCounts, orbitalCounts);
    this->bcast_ScData(sc_file, this->get_nat(), this->get_ntype());
    this->set_npol(NPOL);
    this->set_ParaV(ParaV_in);
    this->set_solver_parameters(nspin_in, kv_in, phsol_in, p_hamilt_in, psi_in, pelec_in, KS_SOLVER_in, LM_in);
}