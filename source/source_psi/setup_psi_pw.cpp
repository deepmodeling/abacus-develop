#include "source_psi/setup_psi_pw.h"
#include "source_io/module_parameter/parameter.h" // use parameter

Setup_Psi_pw::Setup_Psi_pw(){}

Setup_Psi_pw::~Setup_Psi_pw(){}

template <typename T, typename Device>
void Setup_Psi_pw::before_runner_impl(
        const UnitCell &ucell,
        const K_Vectors &kv,
        const Structure_Factor &sf,
        const ModulePW::PW_Basis_K &pw_wfc, 
        const pseudopot_cell_vnl &ppcell,
        const Input_para &inp)
{
    this->p_psi_init = new psi::PSIPrepare<T, Device>(inp.init_wfc,
      inp.ks_solver, inp.basis_type, GlobalV::MY_RANK, ucell,
      sf, kv, ppcell, pw_wfc);

    allocate_psi(this->psi_cpu, kv.get_nks(), kv.ngk, PARAM.globalv.nbands_l, pw_wfc.npwk_max);

    auto* p_psi_init = static_cast<psi::PSIPrepare<T, Device>*>(this->p_psi_init);
    p_psi_init->prepare_init(inp.pw_seed);

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

    if (inp.device == "gpu" || inp.precision == "single") {
        this->psi_t = static_cast<void*>(new psi::Psi<T, Device>(this->psi_cpu[0]));
    } else {
        this->psi_t = static_cast<void*>(reinterpret_cast<psi::Psi<T, Device>*>(this->psi_cpu));
    }
}

void Setup_Psi_pw::before_runner(
        const UnitCell &ucell,
        const K_Vectors &kv,
        const Structure_Factor &sf,
        const ModulePW::PW_Basis_K &pw_wfc, 
        const pseudopot_cell_vnl &ppcell,
        const Input_para &inp)
{
    const bool is_gpu = (inp.device == "gpu");
    const bool is_single = (inp.precision == "single");

#if ((defined __CUDA) || (defined __ROCM))
    if (is_gpu) {
        if (is_single) {
            before_runner_impl<std::complex<float>, base_device::DEVICE_GPU>(
                ucell, kv, sf, pw_wfc, ppcell, inp);
        } else {
            before_runner_impl<std::complex<double>, base_device::DEVICE_GPU>(
                ucell, kv, sf, pw_wfc, ppcell, inp);
        }
    } else
#endif
    {
        if (is_single) {
            before_runner_impl<std::complex<float>, base_device::DEVICE_CPU>(
                ucell, kv, sf, pw_wfc, ppcell, inp);
        } else {
            before_runner_impl<std::complex<double>, base_device::DEVICE_CPU>(
                ucell, kv, sf, pw_wfc, ppcell, inp);
        }
    }
}


template <typename T, typename Device>
void Setup_Psi_pw::update_psi_d_impl()
{
    if (this->psi_d != nullptr && this->precision_type_ == PrecisionType::ComplexFloat)
    {
        delete this->get_psi_d<T, Device>();
    }

    // Refresh this->psi_d
    if (this->precision_type_ == PrecisionType::ComplexFloat) {
        this->psi_d = static_cast<void*>(new psi::Psi<std::complex<double>, Device>(*this->get_psi_t<T, Device>()));
    } else {
        this->psi_d = static_cast<void*>(reinterpret_cast<psi::Psi<std::complex<double>, Device>*>(this->psi_t));
    }
}

void Setup_Psi_pw::update_psi_d()
{
#if ((defined __CUDA) || (defined __ROCM))
    if (this->device_type_ == base_device::GpuDevice)
    {
        if (this->precision_type_ == PrecisionType::ComplexFloat)
        {
            update_psi_d_impl<std::complex<float>, base_device::DEVICE_GPU>();
        }
        else
        {
            update_psi_d_impl<std::complex<double>, base_device::DEVICE_GPU>();
        }
    }
    else
#endif
    {
        if (this->precision_type_ == PrecisionType::ComplexFloat)
        {
            update_psi_d_impl<std::complex<float>, base_device::DEVICE_CPU>();
        }
        else
        {
            update_psi_d_impl<std::complex<double>, base_device::DEVICE_CPU>();
        }
    }
}

template <typename T, typename Device>
void Setup_Psi_pw::init_impl(hamilt::Hamilt<T, Device>* p_hamilt)
{
    if (!this->already_initpsi)
    {
        auto* p_psi_init = static_cast<psi::PSIPrepare<T, Device>*>(this->p_psi_init);
        p_psi_init->initialize_psi(this->psi_cpu, this->get_psi_t<T, Device>(), p_hamilt, GlobalV::ofs_running);
        this->already_initpsi = true;
    }
}

void Setup_Psi_pw::init(hamilt::HamiltBase* p_hamilt)
{
    if (this->already_initpsi)
    {
        return;
    }

#if ((defined __CUDA) || (defined __ROCM))
    if (this->device_type_ == base_device::GpuDevice)
    {
        if (this->precision_type_ == PrecisionType::ComplexFloat)
        {
            init_impl<std::complex<float>, base_device::DEVICE_GPU>(
                static_cast<hamilt::Hamilt<std::complex<float>, base_device::DEVICE_GPU>*>(p_hamilt));
        }
        else
        {
            init_impl<std::complex<double>, base_device::DEVICE_GPU>(
                static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(p_hamilt));
        }
    }
    else
#endif
    {
        if (this->precision_type_ == PrecisionType::ComplexFloat)
        {
            init_impl<std::complex<float>, base_device::DEVICE_CPU>(
                static_cast<hamilt::Hamilt<std::complex<float>, base_device::DEVICE_CPU>*>(p_hamilt));
        }
        else
        {
            init_impl<std::complex<double>, base_device::DEVICE_CPU>(
                static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(p_hamilt));
        }
    }
}


// Transfer data from GPU to CPU in pw basis
template <typename T, typename Device>
void Setup_Psi_pw::copy_d2h_impl()
{
    auto* psi_t = this->get_psi_t<T, Device>();
    this->castmem_d2h_impl<T, Device>(this->psi_cpu[0].get_pointer() - this->psi_cpu[0].get_psi_bias(),
                                      psi_t->get_pointer() - psi_t->get_psi_bias(),
                                      this->psi_cpu[0].size());
}

void Setup_Psi_pw::copy_d2h()
{
    if (this->device_type_ != base_device::GpuDevice)
    {
        return;
    }

#if ((defined __CUDA) || (defined __ROCM))
    if (this->precision_type_ == PrecisionType::ComplexFloat)
    {
        copy_d2h_impl<std::complex<float>, base_device::DEVICE_GPU>();
    }
    else
    {
        copy_d2h_impl<std::complex<double>, base_device::DEVICE_GPU>();
    }
#endif
}

template <typename T, typename Device>
void Setup_Psi_pw::castmem_d2h_impl(std::complex<double>* dst, const std::complex<double>* src, const size_t size)
{
    base_device::memory::cast_memory_op<std::complex<double>, std::complex<double>, base_device::DEVICE_CPU, Device>()(dst, src, size);
}

template <typename T, typename Device>
void Setup_Psi_pw::castmem_d2h_impl(std::complex<double>* dst, const std::complex<float>* src, const size_t size)
{
    base_device::memory::cast_memory_op<std::complex<double>, std::complex<float>, base_device::DEVICE_CPU, Device>()(dst, src, size);
}

template <typename T, typename Device>
void Setup_Psi_pw::clean_impl()
{
    if (this->device_type_ == base_device::GpuDevice || this->precision_type_ == PrecisionType::ComplexFloat)
    {
        delete this->get_psi_t<T, Device>();
    }
    if (this->precision_type_ == PrecisionType::ComplexFloat)
    {
        delete this->get_psi_d<T, Device>();
    }

    delete this->psi_cpu;
    delete this->p_psi_init;
}

void Setup_Psi_pw::clean()
{
#if ((defined __CUDA) || (defined __ROCM))
    if (this->device_type_ == base_device::GpuDevice)
    {
        if (this->precision_type_ == PrecisionType::ComplexFloat)
        {
            clean_impl<std::complex<float>, base_device::DEVICE_GPU>();
        }
        else
        {
            clean_impl<std::complex<double>, base_device::DEVICE_GPU>();
        }
    }
    else
#endif
    {
        if (this->precision_type_ == PrecisionType::ComplexFloat)
        {
            clean_impl<std::complex<float>, base_device::DEVICE_CPU>();
        }
        else
        {
            clean_impl<std::complex<double>, base_device::DEVICE_CPU>();
        }
    }
}

template class psi::PSIPrepare<std::complex<float>, base_device::DEVICE_CPU>;
template class psi::PSIPrepare<std::complex<double>, base_device::DEVICE_CPU>;

template void Setup_Psi_pw::before_runner_impl<std::complex<float>, base_device::DEVICE_CPU>(
    const UnitCell&, const K_Vectors&, const Structure_Factor&,
    const ModulePW::PW_Basis_K&, const pseudopot_cell_vnl&, const Input_para&);

template void Setup_Psi_pw::before_runner_impl<std::complex<double>, base_device::DEVICE_CPU>(
    const UnitCell&, const K_Vectors&, const Structure_Factor&,
    const ModulePW::PW_Basis_K&, const pseudopot_cell_vnl&, const Input_para&);

template void Setup_Psi_pw::init_impl<std::complex<float>, base_device::DEVICE_CPU>(
    hamilt::Hamilt<std::complex<float>, base_device::DEVICE_CPU>*);

template void Setup_Psi_pw::init_impl<std::complex<double>, base_device::DEVICE_CPU>(
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*);

template void Setup_Psi_pw::update_psi_d_impl<std::complex<float>, base_device::DEVICE_CPU>();

template void Setup_Psi_pw::update_psi_d_impl<std::complex<double>, base_device::DEVICE_CPU>();

template void Setup_Psi_pw::clean_impl<std::complex<float>, base_device::DEVICE_CPU>();

template void Setup_Psi_pw::clean_impl<std::complex<double>, base_device::DEVICE_CPU>();

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<float>, base_device::DEVICE_CPU>(
    std::complex<double>*, const std::complex<float>*, const size_t);

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<float>, base_device::DEVICE_CPU>(
    std::complex<double>*, const std::complex<double>*, const size_t);

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<double>, base_device::DEVICE_CPU>(
    std::complex<double>*, const std::complex<float>*, const size_t);

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<double>, base_device::DEVICE_CPU>(
    std::complex<double>*, const std::complex<double>*, const size_t);

#if ((defined __CUDA) || (defined __ROCM))
template class psi::PSIPrepare<std::complex<float>, base_device::DEVICE_GPU>;
template class psi::PSIPrepare<std::complex<double>, base_device::DEVICE_GPU>;

template void Setup_Psi_pw::before_runner_impl<std::complex<float>, base_device::DEVICE_GPU>(
    const UnitCell&, const K_Vectors&, const Structure_Factor&,
    const ModulePW::PW_Basis_K&, const pseudopot_cell_vnl&, const Input_para&);

template void Setup_Psi_pw::before_runner_impl<std::complex<double>, base_device::DEVICE_GPU>(
    const UnitCell&, const K_Vectors&, const Structure_Factor&,
    const ModulePW::PW_Basis_K&, const pseudopot_cell_vnl&, const Input_para&);

template void Setup_Psi_pw::init_impl<std::complex<float>, base_device::DEVICE_GPU>(
    hamilt::Hamilt<std::complex<float>, base_device::DEVICE_GPU>*);

template void Setup_Psi_pw::init_impl<std::complex<double>, base_device::DEVICE_GPU>(
    hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*);

template void Setup_Psi_pw::update_psi_d_impl<std::complex<float>, base_device::DEVICE_GPU>();

template void Setup_Psi_pw::update_psi_d_impl<std::complex<double>, base_device::DEVICE_GPU>();

template void Setup_Psi_pw::copy_d2h_impl<std::complex<float>, base_device::DEVICE_GPU>();

template void Setup_Psi_pw::copy_d2h_impl<std::complex<double>, base_device::DEVICE_GPU>();

template void Setup_Psi_pw::clean_impl<std::complex<float>, base_device::DEVICE_GPU>();

template void Setup_Psi_pw::clean_impl<std::complex<double>, base_device::DEVICE_GPU>();

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<float>, base_device::DEVICE_GPU>(
    std::complex<double>*, const std::complex<float>*, const size_t);

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<float>, base_device::DEVICE_GPU>(
    std::complex<double>*, const std::complex<double>*, const size_t);

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<double>, base_device::DEVICE_GPU>(
    std::complex<double>*, const std::complex<float>*, const size_t);

template void Setup_Psi_pw::castmem_d2h_impl<std::complex<double>, base_device::DEVICE_GPU>(
    std::complex<double>*, const std::complex<double>*, const size_t);
#endif
