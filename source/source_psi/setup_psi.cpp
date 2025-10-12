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

    //! Allocate memory for cpu version of psi
    allocate_psi(this->psi_cpu, this->kv.get_nks(), this->kv.ngk, PARAM.globalv.nbands_l, this->pw_wfc->npwk_max);

    this->p_psi_init->prepare_init(inp.pw_seed);

    //! If GPU or single precision, allocate a new psi (psi_t).
    //! otherwise, transform psi_cpu to psi_t
    this->psi_t = inp.device == "gpu" || inp.precision == "single"
                         ? new psi::Psi<T, Device>(this->psi_cpu[0])
                         : reinterpret_cast<psi::Psi<T, Device>*>(this->psi_cpu);
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
        this->p_psi_init->initialize_psi(this->psi_cpu, this->psi_t, p_hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }
}


template <typename T, typename Device>
void Setup_Psi<T>::copy_g22()
{
    // Transfer data from GPU to CPU in pw basis
    if (this->device == base_device::GpuDevice)
    {
        castmem_2d_d2h_op()(this->psi_cpu[0].get_pointer() - this->psi_cpu[0].get_psi_bias(),
                            this->psi_t[0].get_pointer() - this->psi_t[0].get_psi_bias(),
                            this->psi_cpu[0].size());
    }
}



template <typename T, typename Device>
void Setup_Psi<T>::clean()
{
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        delete this->psi_t;
    }
    if (PARAM.inp.precision == "single")
    {
        delete this->psi_d;
    }

    delete this->psi_cpu;
    delete this->p_psi_init;
}
