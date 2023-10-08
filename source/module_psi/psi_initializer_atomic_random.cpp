#include "psi_initializer_atomic_random.h"

#ifdef __MPI
psi_initializer_atomic_random::psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in) 
                              : psi_initializer_atomic(sf_in, pw_wfc_in, p_ucell_in, p_parakpts_in, random_seed_in)
#else
psi_initializer_atomic_random::psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in) 
                              : psi_initializer_atomic(sf_in, pw_wfc_in, p_ucell_in, random_seed_in)
#endif
{
    this->set_method("atomic+random");
    this->set_random_mix(0.05);
}

psi_initializer_atomic_random::~psi_initializer_atomic_random() {}

void psi_initializer_atomic_random::initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in)
{
    psi_initializer_atomic::initialize_only_once(p_pspot_nl_in);
}

psi::Psi<std::complex<double>>* psi_initializer_atomic_random::cal_psig(int ik)
{
    double rm = this->get_random_mix();
    this->psig->fix_k(ik);
    this->psig = psi_initializer_atomic::cal_psig(ik);
    psi::Psi<std::complex<double>> psi_random(1, this->psig->get_nbands(), this->psig->get_nbasis(), nullptr);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, psi_random.get_nbands(), ik);
    for(int iband = 0; iband < this->psig->get_nbands(); iband++)
    {
        for(int ibasis = 0; ibasis < this->psig->get_nbasis(); ibasis++)
        {
            (*(this->psig))(iband, ibasis) = (1-rm)*(*(this->psig))(iband, ibasis) + rm*psi_random(iband, ibasis);
        }
    }
    return this->psig;
}