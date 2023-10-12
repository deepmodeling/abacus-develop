// representation transform
#include "module_basis/module_representation/representation.h"
// input representation
#include "module_basis/module_representation/repin_pw.h"       // for pw restart
#include "module_basis/module_representation/repin_pao.h"      // for init_wfc atomic, atomic+random
#include "module_basis/module_representation/repin_nao.h"      // for init_wfc nao, nao+random
// output representation
#include "module_basis/module_representation/repout_pw.h"      // for COHP
#include "module_basis/module_representation/repout_qo.h"      // for DeepTB
#include "module_basis/module_representation/repout_wannier.h" // for DeepTB
// for RepIn and RepOut initialization
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
//#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
//#endif
template<typename T, typename Device>
Representation<T, Device>::Representation()
{
    this->repin = nullptr;
    this->representations["input"] = {};
    this->representations["output"] = {};
}

template<typename T, typename Device>
Representation<T, Device>::~Representation()
{
    if (this->repin != nullptr)
    {
        delete this->repin;
    }
    for (auto& r : this->repout)
    {
        if (r != nullptr) // but actually this judgement is not necessary, in principle there is no possibility to push nullptr into it
        {
            delete r;
        }
    }
}

template<typename T, typename Device>
void Representation<T, Device>::transform(const psi::Psi<T, Device>& psi_in, 
                                                psi::Psi<T, Device>& psi_out)
{
    if (this->repin == nullptr)
    {
        throw std::runtime_error("Representation::transform: repin is not set");
    }
    this->repin->cal_psig(this->psig);
    for (auto& r : this->repout)
    {
        //        nbands*nbasis1 nbasis1*nbasis2 nbands*nbasis2
        // LCAO:  nbands*nlocal  nlocal*npwx     nbands*npwx
        r->project(psi_in,       this->psig,     psi_out);
    }
}

template<typename T, typename Device>
void Representation<T, Device>::configure(Structure_Factor* sf_in, 
                                          ModulePW::PW_Basis_K* pw_wfc_in, 
                                          UnitCell* p_ucell_in, 
                                          Parallel_Kpoints* p_parakpts_in, 
                                          pseudopot_cell_vnl* p_pspot_nl_in)
{
    this->p_sf = sf_in;
    this->pw_wfc = pw_wfc_in;
    this->p_ucell = p_ucell_in;
    this->p_parakpts = p_parakpts_in;
    this->p_pspot_nl = p_pspot_nl_in;
}

template<typename T, typename Device>
void Representation<T, Device>::set_repin(const std::string repin_name)
{
    if (this->repin != nullptr)
    {
        delete this->repin;
    }
    if (repin_name == "pw")
    {
        this->repin = new RepIn_PW<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_parakpts, this->p_pspot_nl);
    }
    else if (repin_name == "pao")
    {
        this->repin = new RepIn_PAO<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_parakpts, this->p_pspot_nl);
    }
    else if (repin_name == "nao")
    {
        this->repin = new RepIn_NAO<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_parakpts, this->p_pspot_nl);
    }
    else
    {
        throw std::runtime_error("Representation<T, Device>::set_repin: repin_name is not supported");
    }
    this->representations["input"][0] = repin_name;
}

template<typename T, typename Device>
void Representation<T, Device>::add_repout(const std::string repout_name)
{
    if (repout_name == "pw")
    {
        this->repout.push_back(new RepOut_PW<T, Device>());
    }
    else if (repout_name == "qo")
    {
        this->repout.push_back(new RepOut_QO<T, Device>());
    }
    else if (repout_name == "wannier")
    {
        this->repout.push_back(new RepOut_Wannier<T, Device>());
    }
    else
    {
        throw std::runtime_error("Representation<T, Device>::add_repout: repout_name is not supported");
    }
    this->repout.back()->set_rep_from(this->representations["input"][0]);
    this->representations["output"].push_back(repout_name);
}

template<typename T, typename Device>
void Representation<T, Device>::add_transform_pair(const std::string repin_name, 
                                                   const std::string repout_name)
{
    if(this->representations["input"].size() != 0)
    {
        if(this->representations["input"][0] != repin_name)
        {
            throw std::runtime_error("Representation<T, Device>::add_transform_pair: repin_name is not consistent");
        }
        else
        {
            // do nothing, because repin is already set and allocated
        }
    }
    else
    {
        this->set_repin(repin_name);
    }
    this->add_repout(repout_name);
}

template<typename T, typename Device>
void Representation<T, Device>::clean_representations()
{
    // clean representation name
    this->representations["input"].clear();
    this->representations["output"].clear();

    // deallocate repin
    if (this->repin != nullptr)
    {
        delete this->repin;
        this->repin = nullptr;
    }
    // deallocate all repout
    for (auto& r : this->repout)
    {
        if (r != nullptr)
        {
            delete r;
            r = nullptr;
        }
    }

}

// explicit instantiation
template class Representation<std::complex<double>, psi::DEVICE_CPU>;
template class Representation<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation support
template class Representation<double, psi::DEVICE_CPU>;
template class Representation<float, psi::DEVICE_CPU>;
// GPU support
#if ((defined __CUDA) || (defined __ROCM))
template class Representation<std::complex<double>, psi::DEVICE_GPU>;
template class Representation<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation support
template class Representation<double, psi::DEVICE_GPU>;
template class Representation<float, psi::DEVICE_GPU>;
#endif