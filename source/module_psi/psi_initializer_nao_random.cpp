#include "psi_initializer_nao_random.h"

#ifdef __MPI
template <typename T, typename Device>
void psi_initializer_nao_random<T, Device>::initialize(Structure_Factor* sf,
                                                       ModulePW::PW_Basis_K* pw_wfc,
                                                       UnitCell* p_ucell,
                                                       Parallel_Kpoints* p_parakpts,
                                                       const int& random_seed,
                                                       pseudopot_cell_vnl* p_pspot_nl,
                                                       const int& rank)
{
    psi_initializer_nao<T, Device>::initialize(sf, pw_wfc, p_ucell, p_parakpts, random_seed, p_pspot_nl, rank);
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
}
#else
template <typename T, typename Device>
void psi_initializer_nao_random<T, Device>::initialize(Structure_Factor* sf,
                                                       ModulePW::PW_Basis_K* pw_wfc,
                                                       UnitCell* p_ucell,
                                                       const int& random_seed,
                                                       pseudopot_cell_vnl* p_pspot_nl)
{
    psi_initializer_nao<T, Device>::initialize(sf, pw_wfc, p_ucell, random_seed, p_pspot_nl);
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
}
#endif

template <typename T, typename Device>
psi::Psi<T, Device>* psi_initializer_nao_random<T, Device>::cal_psig(int ik)
{
    double rm = this->random_mix();
    this->psig->fix_k(ik);
    this->psig = psi_initializer_nao<T, Device>::cal_psig(ik);
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

template class psi_initializer_nao_random<std::complex<double>, psi::DEVICE_CPU>;
template class psi_initializer_nao_random<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer_nao_random<double, psi::DEVICE_CPU>;
template class psi_initializer_nao_random<float, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer_nao_random<std::complex<double>, psi::DEVICE_GPU>;
template class psi_initializer_nao_random<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer_nao_random<double, psi::DEVICE_GPU>;
template class psi_initializer_nao_random<float, psi::DEVICE_GPU>;
#endif