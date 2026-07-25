#include "psi.h"
#include "source_base/tool_quit.h"

#include <cstring>
#include <type_traits>

namespace psi
{

namespace detail
{
template <typename T, typename Device>
struct PagingTransfer
{
    static void load(T* dst, const T* src, size_t size);
    static void store(T* dst, const T* src, size_t size);
};

template <typename T>
struct PagingTransfer<T, base_device::DEVICE_CPU>
{
    static void load(T* dst, const T* src, size_t size)
    {
        std::memcpy(dst, src, sizeof(T) * size);
    }
    static void store(T* dst, const T* src, size_t size)
    {
        std::memcpy(dst, src, sizeof(T) * size);
    }
};

#if defined(__CUDA) || defined(__ROCM)
template <typename T>
struct PagingTransfer<T, base_device::DEVICE_GPU>
{
    static void load(T* dst, const T* src, size_t size)
    {
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(dst, src, size);
    }
    static void store(T* dst, const T* src, size_t size)
    {
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(dst, src, size);
    }
};
#endif
}

template <typename T, typename Device>
void Psi<T, Device>::set_storage_mode(PsiStorageMode mode)
{
    if (mode == storage_mode_)
    {
        return;
    }

    if (this->psi != nullptr)
    {
        if (storage_mode_ == PsiStorageMode::ALL_GPU && mode == PsiStorageMode::PAGED_GPU)
        {
            ModuleBase::WARNING("Psi::set_storage_mode",
                                "Setting storage mode flag on already-allocated Psi for copy construction propagation.");
        }
        else
        {
            ModuleBase::WARNING_QUIT("Psi::set_storage_mode",
                                     "Cannot change storage mode after allocation");
        }
    }

    storage_mode_ = mode;
}

template <typename T, typename Device>
T* Psi<T, Device>::get_cpu_pointer(int ik)
{
    if (storage_mode_ == PsiStorageMode::PAGED_GPU)
    {
        if (ik < 0 || ik >= this->nk)
        {
            ModuleBase::WARNING_QUIT("Psi::get_cpu_pointer", "Invalid k-point index");
        }
        if (psi_cpu_ == nullptr)
        {
            ModuleBase::WARNING_QUIT("Psi::get_cpu_pointer",
                                     "CPU buffer not allocated in PAGED_GPU mode");
        }
        return psi_cpu_ + static_cast<size_t>(ik) * this->nbands * this->nbasis;
    }
    else
    {
        return this->psi + static_cast<size_t>(ik) * this->nbands * this->nbasis;
    }
}

template <typename T, typename Device>
const T* Psi<T, Device>::get_cpu_pointer(int ik) const
{
    if (storage_mode_ == PsiStorageMode::PAGED_GPU)
    {
        if (ik < 0 || ik >= this->nk)
        {
            ModuleBase::WARNING_QUIT("Psi::get_cpu_pointer", "Invalid k-point index");
        }
        if (psi_cpu_ == nullptr)
        {
            ModuleBase::WARNING_QUIT("Psi::get_cpu_pointer",
                                     "CPU buffer not allocated in PAGED_GPU mode");
        }
        return psi_cpu_ + static_cast<size_t>(ik) * this->nbands * this->nbasis;
    }
    else
    {
        return this->psi + static_cast<size_t>(ik) * this->nbands * this->nbasis;
    }
}

template <typename T, typename Device>
const T* Psi<T, Device>::get_cpu_pointer_safe(int ik) const
{
    if (storage_mode_ == PsiStorageMode::PAGED_GPU && psi_cpu_ != nullptr)
    {
        if (ik < 0 || ik >= this->nk) return nullptr;
        return psi_cpu_ + static_cast<size_t>(ik) * this->nbands * this->nbasis;
    }
    return nullptr;
}

template <typename T, typename Device>
void Psi<T, Device>::load_k_to_gpu(int ik)
{
    if (storage_mode_ != PsiStorageMode::PAGED_GPU)
    {
        return;
    }

    if (ik < 0 || ik >= this->nk)
    {
        ModuleBase::WARNING_QUIT("Psi::load_k_to_gpu", "Invalid k-point index");
    }

    if (psi_cpu_ == nullptr)
    {
        ModuleBase::WARNING_QUIT("Psi::load_k_to_gpu", "CPU buffer not allocated");
    }

    const size_t k_size = static_cast<size_t>(this->nbands) * this->nbasis;
    const T* src = psi_cpu_ + static_cast<size_t>(ik) * this->nbands * this->nbasis;

    detail::PagingTransfer<T, Device>::load(this->psi, src, k_size);

    current_k_gpu_ = ik;
    this->psi_current = this->psi;
}

template <typename T, typename Device>
void Psi<T, Device>::store_k_from_gpu(int ik)
{
    if (storage_mode_ != PsiStorageMode::PAGED_GPU)
    {
        return;
    }

    if (ik < 0 || ik >= this->nk)
    {
        ModuleBase::WARNING_QUIT("Psi::store_k_from_gpu", "Invalid k-point index");
    }

    if (psi_cpu_ == nullptr)
    {
        ModuleBase::WARNING_QUIT("Psi::store_k_from_gpu", "CPU buffer not allocated");
    }

    const size_t k_size = static_cast<size_t>(this->nbands) * this->nbasis;
    T* dst = psi_cpu_ + static_cast<size_t>(ik) * this->nbands * this->nbasis;

    detail::PagingTransfer<T, Device>::store(dst, this->psi, k_size);
}

template <typename T, typename Device>
void Psi<T, Device>::ensure_k_on_gpu(int ik)
{
    if (storage_mode_ != PsiStorageMode::PAGED_GPU)
    {
        return;
    }

    if (current_k_gpu_ != ik)
    {
        load_k_to_gpu(ik);
    }
}

template <typename T, typename Device>
void Psi<T, Device>::set_psi_cpu_external(T* ext_cpu_buf)
{
    if (storage_mode_ != PsiStorageMode::PAGED_GPU)
    {
        ModuleBase::WARNING_QUIT("Psi::set_psi_cpu_external",
                                 "Only valid in PAGED_GPU mode");
    }
    if (psi_cpu_ != nullptr && psi_cpu_owned_)
    {
        delete[] psi_cpu_;
    }
    psi_cpu_ = ext_cpu_buf;
    psi_cpu_owned_ = false;
}

#define INSTANTIATE_PAGING_METHODS(T, Device)                                                                          \
    template void Psi<T, Device>::set_storage_mode(PsiStorageMode);                                                    \
    template T* Psi<T, Device>::get_cpu_pointer(int);                                                                  \
    template const T* Psi<T, Device>::get_cpu_pointer(int) const;                                                      \
    template const T* Psi<T, Device>::get_cpu_pointer_safe(int) const;                                                 \
    template void Psi<T, Device>::load_k_to_gpu(int);                                                                  \
    template void Psi<T, Device>::store_k_from_gpu(int);                                                               \
    template void Psi<T, Device>::ensure_k_on_gpu(int);                                                                \
    template void Psi<T, Device>::set_psi_cpu_external(T*);

INSTANTIATE_PAGING_METHODS(float, base_device::DEVICE_CPU)
INSTANTIATE_PAGING_METHODS(std::complex<float>, base_device::DEVICE_CPU)
INSTANTIATE_PAGING_METHODS(double, base_device::DEVICE_CPU)
INSTANTIATE_PAGING_METHODS(std::complex<double>, base_device::DEVICE_CPU)
#if ((defined __CUDA) || (defined __ROCM))
INSTANTIATE_PAGING_METHODS(float, base_device::DEVICE_GPU)
INSTANTIATE_PAGING_METHODS(std::complex<float>, base_device::DEVICE_GPU)
INSTANTIATE_PAGING_METHODS(double, base_device::DEVICE_GPU)
INSTANTIATE_PAGING_METHODS(std::complex<double>, base_device::DEVICE_GPU)
#endif

#undef INSTANTIATE_PAGING_METHODS

} // namespace psi
