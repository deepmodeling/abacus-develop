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
#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#endif

#include "module_base/timer.h"
#include "module_base/tool_quit.h"

template<typename T, typename Device>
Representation<T, Device>::Representation()
{
    this->representations["input"].clear();
    this->representations["input"].shrink_to_fit();
    this->representations["output"].clear();
    this->representations["output"].shrink_to_fit();
}

template<typename T, typename Device>
Representation<T, Device>::~Representation()
{
    this->clean_representations();
    this->deallocate_psig();
}

template<typename T, typename Device>
void Representation<T, Device>::transform(psi::Psi<T, Device>* psi_in, 
                                          psi::Psi<T, Device>* psi_out)
{
    ModuleBase::timer::tick("Representation", "transform");
    if (this->repin == nullptr)
    {
        ModuleBase::WARNING_QUIT("Representation::transform", "repin is not allocated");
    }
    // it is one of central function in the whole representation module, pw representation of basis functions is calculated here
    // then if input wavefunction is not nullptr, will multiply it with psig to get the transformed wavefunction and write to
    // psi_out.
    this->repin->cal_psig(this->psig);
    // because psi always occupies large amount of memory, so it is not recommended allocate a lot of psi in different representation
    // and then do transform, instead, do transform one by one, and then deallocate psi.
    // For the demand of multiple representation psi output, can design such a loop:
    // (if there are keywords like out_qo, out_cohp, out_wannier, out_gto), then there will be a vector storing these values and,
    // suppose the vector is std::vector<std::string> out_representations
    // ```C++
    //     while(out_representations.size() > 0)
    //     {
    //         this->rep.clean_representations();
    //         this->rep.add_transform_pair(this->representations["input"][0], out_representations[0]);
    //         ...
    //         this->transform(psi_in, psi_out);
    // ```
    if(this->repout.size() > 1)
    {
        ModuleBase::WARNING_QUIT("Representation::transform", "presently only one repout is supported");
    }
    for (auto& r : this->repout)
    {
        if(psig == nullptr)
        {
            ModuleBase::WARNING_QUIT("Representation::transform", "psig is not allocated");
        }
        else
        {
            if(psi_in == nullptr)
            {
                // it is the case psi_in pretends to be identity matrix
                if(psi_out == nullptr)
                {
                    // it is the case psi_out is not needed to calculate,
                    // we define this case is single-calculation of psig, for pw wavefunction initialization
                    return;
                }
                else
                {
                    // in this case we have psi_out, but we dont have psi_in, or say psi_in is just a identity matrix
                    // directly copy elements from psig to psi_out, but care about the size
                    int nbands = psig->get_nbands()>psi_out->get_nbands()?psi_out->get_nbands():psig->get_nbands();
                    int nbasis = psig->get_nbasis()>psi_out->get_nbasis()?psi_out->get_nbasis():psig->get_nbasis();
                    for(int iband = 0; iband < nbands; iband++)
                    {
                        for(int ibasis = 0; ibasis < nbasis; ibasis++)
                        {
                            psi_out->get_pointer()[iband * psi_out->get_nbasis() + ibasis] = psig->get_pointer()[iband * psig->get_nbasis() + ibasis];
                        }
                    }
                }
            }
            else
            {
                // we really have psi_in in this case
                if(psi_out == nullptr)
                {
                    // but we dont have output, in this case, warning and quit
                    ModuleBase::WARNING_QUIT("Representation::transform", "psi_out is nullptr");
                }
                else
                {
                    // we have psi_in and psi_out, so do gemm to calculate psi_out = psi_in * psig
                    //        psi::Psi&      psi::Psi*       psi::Psi&
                    //        nbands*nbasis1 nbasis1*nbasis2 nbands*nbasis2
                    // LCAO:  nbands*nlocal  nlocal*npwx     nbands*npwx
                    r->project(psi_in,       this->psig,     psi_out);
                }
            }
        }
    }
    ModuleBase::timer::tick("Representation", "transform");
}

template<typename T, typename Device>
#ifdef __MPI
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
#else
void Representation<T, Device>::configure(Structure_Factor* sf_in, 
                                          ModulePW::PW_Basis_K* pw_wfc_in, 
                                          UnitCell* p_ucell_in, 
                                          pseudopot_cell_vnl* p_pspot_nl_in)
{
    this->p_sf = sf_in;
    this->pw_wfc = pw_wfc_in;
    this->p_ucell = p_ucell_in;
    this->p_pspot_nl = p_pspot_nl_in;
}
#endif

template<typename T, typename Device>
void Representation<T, Device>::set_repin(const std::string repin_name)
{
    if (this->repin != nullptr)
    {
        delete this->repin;
        this->repin = nullptr;
    }

    if (repin_name == "pao")
    {
        #ifdef __MPI
        this->repin = new RepIn_PAO<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_parakpts, this->p_pspot_nl);
        #else
        this->repin = new RepIn_PAO<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_pspot_nl);
        #endif
        this->external_files = this->p_ucell->pseudo_fn;
    }
    else if (repin_name == "nao")
    {
        #ifdef __MPI
        this->repin = new RepIn_NAO<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_parakpts, this->p_pspot_nl);
        #else
        this->repin = new RepIn_NAO<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_pspot_nl);
        #endif
        this->external_files = this->p_ucell->orbital_fn;
    }
/*
    else if (repin_name == "pw")
    {
        this->repin = new RepIn_PW<T, Device>(
            this->p_sf, this->pw_wfc, this->p_ucell, this->p_parakpts, this->p_pspot_nl);
    }
*/
    else
    {
        throw std::runtime_error("Representation<T, Device>::set_repin: repin_name is not supported");
    }
    this->representations["input"].push_back(repin_name);
}

template<typename T, typename Device>
void Representation<T, Device>::add_repout(const std::string repout_name)
{
    if (repout_name == "pw")
    {
        this->repout.push_back(new RepOut_PW<T, Device>());
    }
/*
    else if (repout_name == "qo")
    {
        ModuleBase::WARNING_QUIT("Representation<T, Device>::add_repout", "repout_name is not supported");
        this->repout.push_back(new RepOut_QO<T, Device>());
    }
*/
/*
    else if (repout_name == "wannier")
    {
        ModuleBase::WARNING_QUIT("Representation<T, Device>::add_repout", "repout_name is not supported");
        this->repout.push_back(new RepOut_Wannier<T, Device>());
    }
*/
    else
    {
        ModuleBase::WARNING_QUIT("Representation<T, Device>::add_repout", "repout_name is not supported");
    }
    this->representations["output"].push_back(repout_name);
}

template<typename T, typename Device>
void Representation<T, Device>::add_transform_pair(const std::string repin_name, 
                                                   const std::string repout_name)
{
    ModuleBase::timer::tick("Representation", "add_transform_pair");
    if(this->representations["input"].size() != 0)
    {
        if(this->representations["input"][0] != repin_name)
        {
            ModuleBase::WARNING_QUIT("Representation::add_transform_pair", 
            "to change psi want to transform, first do Representation::clean_representation() first.");
        }
        else
        {
            // do nothing, because repin is already set and allocated
        }
    }
    else
    {
        this->set_repin(repin_name);
        this->repin->initialize(this->external_files);
    }
    this->add_repout(repout_name);
    ModuleBase::timer::tick("Representation", "add_transform_pair");
}

template<typename T, typename Device>
void Representation<T, Device>::clean_representations()
{
    // clean representation name
    this->representations["input"].clear();
    this->representations["input"].shrink_to_fit();
    this->representations["output"].clear();
    this->representations["output"].shrink_to_fit();

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
    this->repout.clear();
    this->repout.shrink_to_fit();
}

template<typename T, typename Device>
void Representation<T, Device>::allocate_psig(const int nks, const int nbands, const int nbasis, const int* npwk)
{
    this->psig = new psi::Psi<T, Device>(nks, nbands, nbasis, npwk);
}

template<typename T, typename Device>
void Representation<T, Device>::deallocate_psig()
{
    if (this->psig != nullptr)
    {
        delete this->psig;
        this->psig = nullptr;
    }
}

template<typename T, typename Device>
void Representation<T, Device>::align_kpoint(int ik_in)
{
    this->repin->set_kpoint(ik_in);
    this->psig->fix_k(ik_in);
    for (auto& r : this->repout)
    {
        r->set_kpoint(ik_in);
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