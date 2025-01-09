#include "psi_initializer_file.h"

#include "module_base/timer.h"
#include "module_cell/klist.h"
#include "module_io/read_wfc_pw.h"
#include "module_parameter/parameter.h"
template <typename T>
void psi_initializer_file<T>::initialize(const Structure_Factor* sf,
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
    this->nbands_start_ = PARAM.inp.nbands;
    this->nbands_complem_ = 0;
}

template <typename T>
void psi_initializer_file<T>::init_psig(T* psig, const int& ik)
{
    ModuleBase::timer::tick("psi_initializer_file", "init_psig");
    const int npol = PARAM.globalv.npol;
    const int nbasis = this->pw_wfc_->npwk_max * npol;
    const int pre_nk = PARAM.inp.nspin == 2 ? 2 : 1;
    const int nkstot = this->p_parakpts_->nkstot_np * pre_nk;
    ModuleBase::ComplexMatrix wfcatom(this->nbands_start_, nbasis);
    std::stringstream filename;
    int ik_tot = K_Vectors::get_ik_global(ik, nkstot);
    filename << PARAM.globalv.global_readin_dir << "WAVEFUNC" << ik_tot + 1 << ".dat";
    ModuleIO::read_wfc_pw(filename.str(), this->pw_wfc_, ik, nkstot, wfcatom);

    assert(this->nbands_start_ <= wfcatom.nr);
    for (int ib = 0; ib < this->nbands_start_; ib++)
    {
        for (int ig = 0; ig < nbasis; ig++)
        {
            psig[ib * nbasis + ig] = this->template cast_to_T<T>(wfcatom(ib, ig));
        }
    }
    ModuleBase::timer::tick("psi_initializer_file", "init_psig");
}

template class psi_initializer_file<std::complex<double>>;
template class psi_initializer_file<std::complex<float>>;