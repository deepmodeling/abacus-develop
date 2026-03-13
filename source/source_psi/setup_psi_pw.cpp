#include "source_psi/setup_psi_pw.h"
#include "source_io/module_parameter/parameter.h" // use parameter

template <typename T, typename Device>
Setup_Psi_pw<T, Device>::Setup_Psi_pw(){}

template <typename T, typename Device>
Setup_Psi_pw<T, Device>::~Setup_Psi_pw(){}

template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::before_runner(
        const UnitCell &ucell,
        const K_Vectors &kv,
        const Structure_Factor &sf,
        const ModulePW::PW_Basis_K &pw_wfc, 
        const pseudopot_cell_vnl &ppcell,
        const Input_para &inp)
{
    //! Allocate and initialize psi
    this->p_psi_init = new psi::PSIPrepare<T, Device>(inp.init_wfc,
      inp.ks_solver, inp.basis_type, GlobalV::MY_RANK, ucell,
      sf, kv, ppcell, pw_wfc);

    //! Allocate memory for cpu version of psi
    allocate_psi(this->psi_cpu, kv.get_nks(), kv.ngk, PARAM.globalv.nbands_l, pw_wfc.npwk_max);

    auto* p_psi_init = static_cast<psi::PSIPrepare<T, Device>*>(this->p_psi_init);
    p_psi_init->prepare_init(inp.pw_seed);

    //! Set runtime type information
    if (std::is_same<T, float>::value) {
        precision_type_ = PrecisionType::Float;
    } else if (std::is_same<T, double>::value) {
        precision_type_ = PrecisionType::Double;
    } else if (std::is_same<T, std::complex<float>>::value) {
        precision_type_ = PrecisionType::ComplexFloat;
    } else {
        precision_type_ = PrecisionType::ComplexDouble;
    }
    
    if (std::is_same<Device, base_device::DEVICE_GPU>::value) {
        device_type_ = base_device::GpuDevice;
    } else {
        device_type_ = base_device::CpuDevice;
    }

    //! If GPU or single precision, allocate a new psi (psi_t).
    //! otherwise, transform psi_cpu to psi_t
    if (inp.device == "gpu" || inp.precision == "single") {
        this->psi_t = static_cast<void*>(new psi::Psi<T, Device>(this->psi_cpu[0]));
    } else {
        this->psi_t = static_cast<void*>(reinterpret_cast<psi::Psi<T, Device>*>(this->psi_cpu));
    }
}


template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::update_psi_d()
{
    if (this->psi_d != nullptr && PARAM.inp.precision == "single")
    {
        delete this->get_psi_d();
    }

    // Refresh this->psi_d
    if (PARAM.inp.precision == "single") {
        this->psi_d = static_cast<void*>(new psi::Psi<std::complex<double>, Device>(*this->get_psi_t()));
    } else {
        this->psi_d = static_cast<void*>(reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->psi_t));
    }
}

template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::init(hamilt::HamiltBase* p_hamilt)
{
    //! Initialize wave functions
    if (!this->already_initpsi)
    {
        auto* p_psi_init = static_cast<psi::PSIPrepare<T, Device>*>(this->p_psi_init);
        auto* hamilt = static_cast<hamilt::Hamilt<T, Device>*>(p_hamilt);
        p_psi_init->initialize_psi(this->psi_cpu, this->get_psi_t(), hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }
}


// Transfer data from GPU to CPU in pw basis (runtime version)
template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::copy_d2h(const base_device::DeviceContext* ctx)
{
    if (base_device::get_device_type(ctx) == base_device::GpuDevice)
    {
        auto* psi_t = this->get_psi_t();
        this->castmem_d2h_impl(this->psi_cpu[0].get_pointer() - this->psi_cpu[0].get_psi_bias(),
                               psi_t->get_pointer() - psi_t->get_psi_bias(),
                               this->psi_cpu[0].size());
    }
    else
    {
       // do nothing
    }
    return;
}

template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::castmem_d2h_impl(std::complex<double>* dst, const std::complex<double>* src, const size_t size)
{
    base_device::memory::cast_memory_op<std::complex<double>, std::complex<double>, base_device::DEVICE_CPU, Device>()(dst, src, size);
}

template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::castmem_d2h_impl(std::complex<double>* dst, const std::complex<float>* src, const size_t size)
{
    base_device::memory::cast_memory_op<std::complex<double>, std::complex<float>, base_device::DEVICE_CPU, Device>()(dst, src, size);
}



template <typename T, typename Device>
void Setup_Psi_pw<T, Device>::clean()
{
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        delete this->get_psi_t();
    }
    if (PARAM.inp.precision == "single")
    {
        delete this->get_psi_d();
    }

    delete this->psi_cpu;
    delete this->p_psi_init;
}

template class Setup_Psi_pw<std::complex<float>, base_device::DEVICE_CPU>;
template class Setup_Psi_pw<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Setup_Psi_pw<std::complex<float>, base_device::DEVICE_GPU>;
template class Setup_Psi_pw<std::complex<double>, base_device::DEVICE_GPU>;
#endif
