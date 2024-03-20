#include "psi_initializer_atomic_random.h"

#ifdef __MPI
template <typename T, typename Device>
void psi_initializer_atomic_random<T, Device>::initialize(Structure_Factor* sf,                                //< structure factor
                                                          ModulePW::PW_Basis_K* pw_wfc,                        //< planewave basis
                                                          UnitCell* p_ucell,                                   //< unit cell
                                                          Parallel_Kpoints* p_parakpts,                        //< parallel kpoints
                                                          const int& random_seed,                          //< random seed
                                                          pseudopot_cell_vnl* p_pspot_nl,
                                                          const int& rank)
{
    psi_initializer_atomic<T, Device>::initialize(sf, pw_wfc, p_ucell, p_parakpts, random_seed, p_pspot_nl, rank);
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
}
#else
template <typename T, typename Device>
void psi_initializer_atomic_random<T, Device>::initialize(Structure_Factor* sf,                                //< structure factor
                                                          ModulePW::PW_Basis_K* pw_wfc,                        //< planewave basis
                                                          UnitCell* p_ucell,                                   //< unit cell
                                                          const int& random_seed,                          //< random seed
                                                          pseudopot_cell_vnl* p_pspot_nl)
{
    psi_initializer_atomic<T, Device>::initialize(sf, pw_wfc, p_ucell, random_seed, p_pspot_nl);
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
}
#endif

template <typename T, typename Device>
void psi_initializer_atomic_random<T, Device>::proj_ao_onkG(int ik)
{
    double rm = this->random_mix();
    this->psig_->fix_k(ik);
    psi_initializer_atomic<T, Device>::proj_ao_onkG(ik);
    psi::Psi<T, Device> psi_random(1, this->psig_->get_nbands(), this->psig_->get_nbasis(), nullptr);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, psi_random.get_nbands(), ik);
    for(int iband = 0; iband < this->psig_->get_nbands(); iband++)
    {
        for(int ibasis = 0; ibasis < this->psig_->get_nbasis(); ibasis++)
        {
            (*(this->psig_))(iband, ibasis) = ((Real)(1-rm))*(*(this->psig_))(iband, ibasis) + ((Real)rm)*psi_random(iband, ibasis);
        }
    }
}

template class psi_initializer_atomic_random<std::complex<double>, psi::DEVICE_CPU>;
template class psi_initializer_atomic_random<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer_atomic_random<double, psi::DEVICE_CPU>;
template class psi_initializer_atomic_random<float, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer_atomic_random<std::complex<double>, psi::DEVICE_GPU>;
template class psi_initializer_atomic_random<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer_atomic_random<double, psi::DEVICE_GPU>;
template class psi_initializer_atomic_random<float, psi::DEVICE_GPU>;
#endif
