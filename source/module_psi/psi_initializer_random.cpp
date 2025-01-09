#include "psi_initializer_random.h"
#ifdef __MPI
#include <mpi.h>
#endif
#include "module_base/parallel_global.h"
#include "module_base/timer.h"
#include "module_cell/parallel_kpoints.h"
#include "module_parameter/parameter.h"

template <typename T>
void psi_initializer_random<T>::initialize(const Structure_Factor* sf,
                                           const ModulePW::PW_Basis_K* pw_wfc,
                                           const UnitCell* p_ucell,
                                           const Parallel_Kpoints* p_parakpts,
                                           const int& random_seed,
                                           const pseudopot_cell_vnl* p_pspot_nl,
                                           const int& rank)
{
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->p_parakpts_ = p_parakpts;
    this->random_seed_ = random_seed;
    this->p_pspot_nl_ = p_pspot_nl;
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
    this->nbands_start_ = PARAM.inp.nbands;
    this->nbands_complem_ = 0;
}

template <typename T>
void psi_initializer_random<T>::init_psig(T* psig, const int& ik)
{
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    this->random_t(psig, 0, this->nbands_start_, ik);
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
}

template class psi_initializer_random<std::complex<double>>;
template class psi_initializer_random<std::complex<float>>;
// gamma point calculation
template class psi_initializer_random<double>;
template class psi_initializer_random<float>;