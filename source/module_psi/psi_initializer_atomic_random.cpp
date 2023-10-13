#include "psi_initializer_atomic_random.h"

template <typename T, typename Device>
#ifdef __MPI
psi_initializer_atomic_random<T, Device>::psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in) 
                              : psi_initializer_atomic<T, Device>(sf_in, pw_wfc_in, p_ucell_in, p_parakpts_in, random_seed_in)
#else
psi_initializer_atomic_random<T, Device>::psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in) 
                              : psi_initializer_atomic<T, Device>(sf_in, pw_wfc_in, p_ucell_in, random_seed_in)
#endif
{
    this->set_method("atomic+random");
    this->set_random_mix(0.05);
}

template <typename T, typename Device>
psi_initializer_atomic_random<T, Device>::~psi_initializer_atomic_random() {}

template <typename T, typename Device>
void psi_initializer_atomic_random<T, Device>::initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in)
{
    psi_initializer_atomic<T, Device>::initialize_only_once(p_pspot_nl_in);
}

template <typename T, typename Device>
psi::Psi<T, Device>* psi_initializer_atomic_random<T, Device>::cal_psig(int ik)
{
    double rm = this->get_random_mix();
    this->psig->fix_k(ik);
    this->psig = psi_initializer_atomic<T, Device>::cal_psig(ik);
    psi::Psi<T, Device> psi_random(1, this->psig->get_nbands(), this->psig->get_nbasis(), nullptr);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, psi_random.get_nbands(), ik);
    for(int iband = 0; iband < this->psig->get_nbands(); iband++)
    {
        for(int ibasis = 0; ibasis < this->psig->get_nbasis(); ibasis++)
        {
            (*(this->psig))(iband, ibasis) = ((Real)(1-rm))*(*(this->psig))(iband, ibasis) + ((Real)rm)*psi_random(iband, ibasis);
        }
    }
    return this->psig;
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