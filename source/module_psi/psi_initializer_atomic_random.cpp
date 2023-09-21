#include "psi_initializer_atomic_random.h"

psi_initializer_atomic_random::psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in) : psi_initializer_atomic(sf_in, pw_wfc_in)
{
    this->set_method("atomic+random");
    this->set_random_mix(0.05);
}

psi_initializer_atomic_random::~psi_initializer_atomic_random() {}

psi::Psi<std::complex<double>>* psi_initializer_atomic_random::cal_psig(int ik)
{
    double rm = this->get_random_mix();
    this->psig->fix_k(ik);
    this->psig = psi_initializer_atomic::cal_psig(ik);
    psi::Psi<std::complex<double>> psi_random(1, this->psig->get_nbands(), this->psig->get_nbasis(), nullptr);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, psi_random.get_nbands(), ik, this->pw_wfc);
    for(int iband = 0; iband < this->psig->get_nbands(); iband++)
    {
        for(int ibasis = 0; ibasis < this->psig->get_nbasis(); ibasis++)
        {
            (*(this->psig))(iband, ibasis) = (1-rm)*(*(this->psig))(iband, ibasis) + rm*psi_random(iband, ibasis);
        }
    }
#ifdef PSI_INITIALIZER_TEST
	this->write_psig();
#endif
    return this->psig;
}