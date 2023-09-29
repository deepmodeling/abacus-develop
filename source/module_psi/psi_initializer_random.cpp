#include "psi_initializer_random.h"
#ifdef __MPI
#include <mpi.h>
#endif
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
// basic functions support
#include "module_base/timer.h"

#ifdef __MPI
psi_initializer_random::psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in) 
                       : psi_initializer(sf_in, pw_wfc_in, p_ucell_in, p_parakpts_in, random_seed_in)
#else
psi_initializer_random::psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in) 
                       : psi_initializer(sf_in, pw_wfc_in, p_ucell_in, random_seed_in)
#endif
{
    this->set_method("random");
}

psi_initializer_random::~psi_initializer_random() {}

void psi_initializer_random::random(std::complex<double>* psi,
                                    const int iw_start,
                                    const int iw_end,
                                    const int ik)
{
    ModuleBase::timer::tick("psi_initializer_random", "random");
    this->random_t(psi, iw_start, iw_end, ik);
    ModuleBase::timer::tick("psi_initializer_random", "random");
}

psi::Psi<std::complex<double>>* psi_initializer_random::cal_psig(int ik)
{
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    //this->print_status(psi);
    this->psig->fix_k(ik);
    this->random(this->psig->get_pointer(), 0, this->psig->get_nbands(), ik);
    // we still need to diagonalize the obtained psi from hsolver::DiagoIterAssist::diagH_subspace
    // will do it in HSolver function...
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    return this->psig;
}
