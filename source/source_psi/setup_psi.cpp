#include "source_lcao/setup_deepks.h"
#include "source_lcao/LCAO_domain.h"
#include "source_io/module_parameter/parameter.h" // use parameter

template <typename T, typename Device>
Setup_Psi<T, Device>::Setup_Psi(){}

template <typename T, typename Device>
Setup_Psi<T, Device>::~Setup_Psi(){}

template <typename T, typename Device>
Setup_Psi<T, Device>::before_runner()
{
    //! Allocate and initialize psi
    this->p_psi_init = new psi::PSIInit<T, Device>(inp.init_wfc,
      inp.ks_solver, inp.basis_type, GlobalV::MY_RANK, ucell,
      this->sf, this->kv, this->ppcell, *this->pw_wfc);

    allocate_psi(this->psi, this->kv.get_nks(), this->kv.ngk, PARAM.globalv.nbands_l, this->pw_wfc->npwk_max);

    this->p_psi_init->prepare_init(inp.pw_seed);

    this->kspw_psi = inp.device == "gpu" || inp.precision == "single"
                         ? new psi::Psi<T, Device>(this->psi[0])
                         : reinterpret_cast<psi::Psi<T, Device>*>(this->psi);
}


template <typename T, typename Device>
Setup_Psi<T, Device>::update_psi_d()
{
    if (this->psi_d != nullptr && PARAM.inp.precision == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->psi_t);
    }

    // Refresh this->psi_d
    this->psi_d = PARAM.inp.precision == "single"
                           ? new psi::Psi<std::complex<double>, Device>(this->psi_t[0])
                           : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->psi_t);
}

template <typename T, typename Device>
void Setup_Psi<T>::init()
{
    //! Initialize wave functions
    if (!this->already_initpsi)
    {
        this->p_psi_init->initialize_psi(this->psi, this->kspw_psi, this->p_hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }
}


template <typename T, typename Device>
void Setup_Psi<T>::clean()
{
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        delete this->kspw_psi;
    }
    if (PARAM.inp.precision == "single")
    {
        delete this->__kspw_psi;
    }

    delete this->psi;
    delete this->p_psi_init;


}
