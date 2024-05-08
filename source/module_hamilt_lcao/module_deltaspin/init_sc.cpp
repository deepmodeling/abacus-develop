#include "spin_constrain.h"

// init sc
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::init_sc(const UnitCell& ucell,
                                            int NPOL,
                                            Parallel_Orbitals* ParaV_in,
                                            int nspin_in,
                                            K_Vectors kv_in,
                                            std::string KS_SOLVER_in,
                                            LCAO_Matrix* LM_in,
                                            hsolver::HSolver<FPTYPE, Device>* phsol_in,
                                            hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
                                            psi::Psi<FPTYPE>* psi_in,
                                            elecstate::ElecState* pelec_in)
{
    this->set_atomCounts(ucell.get_atomCounts());
    this->set_orbitalCounts(ucell.get_orbitalCounts());
    this->set_lnchiCounts(ucell.get_lnchiCounts());
    this->set_nspin(nspin_in);
    this->bcast_ScData(this->sc_file, this->get_nat(), this->get_ntype());
    this->set_npol(NPOL);
    this->set_ParaV(ParaV_in);
    this->set_solver_parameters(kv_in, phsol_in, p_hamilt_in, psi_in, pelec_in, KS_SOLVER_in, LM_in);
}

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;
template class SpinConstrain<double, psi::DEVICE_CPU>;