#include "psi_initializer_random.h"
#ifdef __MPI
#include <mpi.h>
#endif
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
// basic functions support
#include "module_base/timer.h"

#ifdef __MPI
template <typename T, typename Device>
void psi_initializer_random<T, Device>::initialize(Structure_Factor* sf,
                                                   ModulePW::PW_Basis_K* pw_wfc,
                                                   UnitCell* p_ucell,
                                                   Parallel_Kpoints* p_parakpts,
                                                   const int& random_seed,
                                                   pseudopot_cell_vnl* p_pspot_nl,
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
}
#else
template <typename T, typename Device>
void psi_initializer_random<T, Device>::initialize(Structure_Factor* sf,
                                                  ModulePW::PW_Basis_K* pw_wfc,
                                                  UnitCell* p_ucell,
                                                  const int& random_seed,
                                                  pseudopot_cell_vnl* p_pspot_nl)
{
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->random_seed_ = random_seed;
    this->p_pspot_nl_ = p_pspot_nl;
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
}
#endif

template <typename T, typename Device>
void psi_initializer_random<T, Device>::random(T* psi,
                                               const int iw_start,
                                               const int iw_end,
                                               const int ik)
{
    ModuleBase::timer::tick("psi_initializer_random", "random");
    this->random_t(psi, iw_start, iw_end, ik);
    ModuleBase::timer::tick("psi_initializer_random", "random");
}

template <typename T, typename Device>
psi::Psi<T, Device>* psi_initializer_random<T, Device>::cal_psig(int ik)
{
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    this->psig->fix_k(ik);
    this->random(this->psig->get_pointer(), 0, this->psig->get_nbands(), ik);
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    return this->psig;
}

template class psi_initializer_random<std::complex<double>, psi::DEVICE_CPU>;
template class psi_initializer_random<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer_random<double, psi::DEVICE_CPU>;
template class psi_initializer_random<float, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer_random<std::complex<double>, psi::DEVICE_GPU>;
template class psi_initializer_random<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer_random<double, psi::DEVICE_GPU>;
template class psi_initializer_random<float, psi::DEVICE_GPU>;
#endif