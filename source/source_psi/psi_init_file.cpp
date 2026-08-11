#include "psi_init_file.h"
#include "source_basis/module_pw/pw_basis_k.h"

#include <vector>
#include <cassert>
#include <string>
#include <complex>

#include "source_base/timer.h"
#include "source_io/module_wf/read_wfc_pw.h"
#include "source_io/module_output/filename.h"

template <typename T>
void psi_init_file<T>::initialize(const Structure_Factor* sf,
                                               const ModulePW::PW_Basis_K* pw_wfc,
                                               const UnitCell* p_ucell,
                                               const std::vector<int>& ik2iktot,
                                               const int& random_seed,
                                               const int& rank,
                                               const int& npol,
                                               const int& nbands)
{
    psi_base<T>::initialize(sf, pw_wfc, p_ucell, ik2iktot, random_seed, rank, npol, nbands);
    this->nbands_start_ = nbands;
    this->nbands_complem_ = 0;
}

template <typename T>
void psi_init_file<T>::prepare_params(const int& nspin,
                                      const std::string& global_readin_dir,
                                      const int& rank_in_pool,
                                      const int& nproc_in_pool,
                                      const int& nkstot)
{
    this->nspin_ = nspin;
    this->global_readin_dir_ = global_readin_dir;
    this->rank_in_pool_ = rank_in_pool;
    this->nproc_in_pool_ = nproc_in_pool;
    this->nkstot_ = nkstot;
}

template <typename T>
void psi_init_file<T>::init_psig(T* psig, const int& ik)
{
    ModuleBase::timer::start("psi_init_file", "init_psig");
    const int npol = this->npol_;
    const int nbasis = this->pw_wfc_->npwk_max * npol;
    const int nkstot = this->nkstot_;
    ModuleBase::ComplexMatrix wfcatom(this->nbands_start_, nbasis);
    int ik_tot = this->ik2iktot_[ik];

    // mohan update, this is for plane wave, 2025-05-17
    const int out_type = 2;
    const bool out_app_flag = false;
    const bool gamma_only = false;
    const int istep = -1;

    std::string fn = ModuleIO::filename_output(this->global_readin_dir_,"wf","pw",
			ik,this->ik2iktot_,this->nspin_,nkstot,
			out_type,out_app_flag,gamma_only,istep);

    ModuleIO::read_wfc_pw(fn, this->pw_wfc_, 
			this->rank_in_pool_, this->nproc_in_pool_,
			this->nbands_start_, this->npol_,
			ik, ik_tot, nkstot, wfcatom);

    assert(this->nbands_start_ <= wfcatom.nr);
    for (int ib = 0; ib < this->nbands_start_; ib++)
    {
        for (int ig = 0; ig < nbasis; ig++)
        {
            psig[ib * nbasis + ig] = this->template cast_to_T<T>(wfcatom(ib, ig));
        }
    }
    ModuleBase::timer::end("psi_init_file", "init_psig");
}

template class psi_init_file<std::complex<double>>;
template class psi_init_file<std::complex<float>>;
