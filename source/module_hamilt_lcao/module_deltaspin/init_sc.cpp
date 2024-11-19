#include "spin_constrain.h"

// init sc
template <typename FPTYPE>
void spinconstrain::SpinConstrain<FPTYPE>::init_sc(double sc_thr_in,
                                            int nsc_in,
                                            int nsc_min_in,
                                            double alpha_trial_in,
                                            double sccut_in,
                                            double sc_drop_thr_in,
                                            const UnitCell& ucell,
                                            Parallel_Orbitals* ParaV_in,
                                            int nspin_in,
                                            K_Vectors& kv_in,
                                            std::string KS_SOLVER_in,
                                            void* p_hamilt_in,
                                            void* psi_in,
                                            elecstate::ElecState* pelec_in)
{
    this->set_input_parameters(sc_thr_in, nsc_in, nsc_min_in, alpha_trial_in, sccut_in, sc_drop_thr_in);
    this->set_atomCounts(ucell.get_atom_Counts());
    this->set_orbitalCounts(ucell.get_orbital_Counts());
    this->set_lnchiCounts(ucell.get_lnchi_Counts());
    this->set_nspin(nspin_in);
    this->set_target_mag(ucell.get_target_mag());
    this->lambda_ = ucell.get_lambda();
    this->constrain_ = ucell.get_constrain();
    this->atomLabels_ = ucell.get_atomLabels();
    this->set_decay_grad();
    if(ParaV_in != nullptr) this->set_ParaV(ParaV_in);
    this->set_solver_parameters(kv_in, p_hamilt_in, psi_in, pelec_in, KS_SOLVER_in);
}

template class spinconstrain::SpinConstrain<std::complex<double>>;
template class spinconstrain::SpinConstrain<double>;