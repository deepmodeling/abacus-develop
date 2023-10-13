#include "module_psi/psi_initializer/psi_initializer_atomic.h"

template<typename T, typename Device>
psi_initializer_atomic<T, Device>::psi_initializer_atomic(
    UnitCell* p_ucell_in,
    ModulePW::PW_Basis_K* pw_wfc_in,
    Parallel_Kpoints* p_parakpts_in,
    Representation<T, Device>* p_rep_in
    ): psi_initializer<T, Device>(
        p_ucell_in,
        pw_wfc_in,
        p_parakpts_in,
        p_rep_in)
{
    // this->rep = rep_in;
}

template<typename T, typename Device>
psi_initializer_atomic<T, Device>::~psi_initializer_atomic()
{
}

template<typename T, typename Device>
void psi_initializer_atomic<T, Device>::initialize(const psi::Psi<T, Device>& psi_out)
{
    // default the psi_in to be identity matrix
    psi::Psi<T, Device> identity = psi::Psi<T, Device>(1, psi_out.get_nbands(), psi_out.get_nbands());
    for(int iband = 0; iband < psi_out.get_nbands(); iband++)
    {
        identity(iband, iband) = static_cast<T>(1.0);
    }
    // do transform
    this->p_rep->add_transform_pair("pao", "pw");
    for(int ik = 0; ik < this->pw_wfc->nks; ik++)
    {
        this->p_rep->align_kpoint(ik);
        this->p_rep->transform(identity, psi_out);
        // if not enough number of bands initialized, fill with random number
        int nbands_complem = this->get_nbands_complem();
        if(nbands_complem != 0)
        {
            int iband_start = GlobalV::NBANDS - nbands_complem;
            int iband_end = GlobalV::NBANDS;
            this->random_t(psi_out->get_pointer(), 
                           iband_start, // range, starting band index (included)
                           iband_end,   // range, ending band index (excluded)
                           ik);         // kpoint index
        }
    }
}