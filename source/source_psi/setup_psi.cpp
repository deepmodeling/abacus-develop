#include "source_psi/setup_psi.h"
#include "source_lcao/setup_deepks.h"
#include "source_lcao/LCAO_domain.h"
#include "source_io/module_parameter/parameter.h" // use parameter

template <typename T, typename Device>
Setup_Psi<T, Device>::Setup_Psi(){}

template <typename T, typename Device>
Setup_Psi<T, Device>::~Setup_Psi(){}

template <typename T, typename Device>
void Setup_Psi<T, Device>::before_runner(
		const UnitCell &ucell,
		const K_Vectors &kv,
		const Structure_Factor &sf,
		const ModulePW::PW_Basis_K &pw_wfc, 
		const pseudopot_cell_vnl &ppcell,
		const Input_para &inp)
{
    //! Allocate and initialize psi
    this->p_psi_init = new psi::PSIInit<T, Device>(inp.init_wfc,
      inp.ks_solver, inp.basis_type, GlobalV::MY_RANK, ucell,
      sf, kv, ppcell, pw_wfc);

    //! Allocate memory for cpu version of psi
    allocate_psi(this->psi_cpu, kv.get_nks(), kv.ngk, PARAM.globalv.nbands_l, pw_wfc.npwk_max);

    this->p_psi_init->prepare_init(inp.pw_seed);

    //! If GPU or single precision, allocate a new psi (psi_t).
    //! otherwise, transform psi_cpu to psi_t
    this->psi_t = inp.device == "gpu" || inp.precision == "single"
                         ? new psi::Psi<T, Device>(this->psi_cpu[0])
                         : reinterpret_cast<psi::Psi<T, Device>*>(this->psi_cpu);
}


template <typename T, typename Device>
void Setup_Psi<T, Device>::update_psi_d()
{
    if (this->psi_d != nullptr && PARAM.inp.precision == "single")
    {
        delete reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->psi_d);
    }

    // Refresh this->psi_d
    this->psi_d = PARAM.inp.precision == "single"
                           ? new psi::Psi<std::complex<double>, Device>(this->psi_t[0])
                           : reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->psi_t);
}

template <typename T, typename Device>
void Setup_Psi<T, Device>::init(hamilt::Hamilt<T, Device>* p_hamilt)
{
    //! Initialize wave functions
    if (!this->already_initpsi)
    {
        this->p_psi_init->initialize_psi(this->psi_cpu, this->psi_t, p_hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }
}


// Transfer data from GPU to CPU in pw basis
template <typename T, typename Device>
void Setup_Psi<T, Device>::copy_d2h(const base_device::AbacusDevice_t &device)
{
    if (device == base_device::GpuDevice)
    {
        castmem_2d_d2h_op()(this->psi_cpu[0].get_pointer() - this->psi_cpu[0].get_psi_bias(),
                            this->psi_t[0].get_pointer() - this->psi_t[0].get_psi_bias(),
                            this->psi_cpu[0].size());
    }
	else
	{
       // do nothing
	}
    return;
}



template <typename T, typename Device>
void Setup_Psi<T, Device>::clean()
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

template class Setup_Psi<std::complex<float>, base_device::DEVICE_CPU>;
template class Setup_Psi<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Setup_Psi<std::complex<float>, base_device::DEVICE_GPU>;
template class Setup_Psi<std::complex<double>, base_device::DEVICE_GPU>;
#endif
