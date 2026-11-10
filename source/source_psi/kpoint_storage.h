#ifndef PSI_KPOINT_STORAGE_H
#define PSI_KPOINT_STORAGE_H

#include "psi.h"

#include <cassert>
#include <stdexcept>
#include <type_traits>
#if defined(__CUDA)
#include <cuda_runtime.h>
namespace psi { namespace detail {
inline void kpoint_check(cudaError_t status)
{
    if (status != cudaSuccess) { throw std::runtime_error(cudaGetErrorString(status)); }
}
} }
#endif

namespace psi
{

/// Temporarily retain only the current k point on Device. Only owned, k-first
/// arrays with all local k points are eligible. The complete device layout is
/// restored after solver workspaces are released, before density or output use.
template <typename T, typename Device>
class KPointStorage
{
  public:
    KPointStorage(Psi<T, Device>& psi, int local_points, bool enabled)
        : psi_(psi), points_(psi.get_nk()), bands_(psi.get_nbands()),
          basis_(psi.get_nbasis()), block_size_(size_t(bands_) * basis_)
    {
#if defined(__CUDA)
        if (!enabled || !std::is_same<Device, base_device::DEVICE_GPU>::value
            || !psi_.allocate_inside || !psi_.k_first || points_ <= 1
            || points_ != local_points || bands_ <= 0 || basis_ <= 0)
        {
            return;
        }
        const size_t bytes = size_t(points_) * block_size_ * sizeof(T);
        const auto status = cudaHostAlloc(reinterpret_cast<void**>(&host_), bytes, cudaHostAllocDefault);
        if (status == cudaErrorMemoryAllocation)
        {
            cudaGetLastError();
            host_ = nullptr;
            return;
        }
        detail::kpoint_check(status);
        detail::kpoint_check(cudaMemcpy(host_, psi_.psi, bytes, cudaMemcpyDeviceToHost));
        const int point = psi_.get_current_k();
        const int band = psi_.get_current_b();
        psi_.resize(1, bands_, basis_);
        psi_.current_k = point;
        psi_.current_b = band;
        psi_.current_nbasis = psi_.ngk != nullptr ? psi_.ngk[point] : basis_;
        psi_.psi_bias = size_t(band) * basis_;
        psi_.psi_current = psi_.psi + psi_.psi_bias;
#endif
    }

    KPointStorage(const KPointStorage&) = delete;
    KPointStorage& operator=(const KPointStorage&) = delete;
    ~KPointStorage() { restore(); }

    bool active() const { return host_ != nullptr; }

    void load(int point)
    {
#if defined(__CUDA)
        if (!active()) { return; }
        assert(point >= 0 && point < points_);
        detail::kpoint_check(cudaMemcpy(psi_.psi, host_ + size_t(point) * block_size_,
                               block_size_ * sizeof(T), cudaMemcpyHostToDevice));
        psi_.current_k = point;
        psi_.current_b = 0;
        psi_.current_nbasis = psi_.ngk != nullptr ? psi_.ngk[point] : basis_;
        psi_.psi_bias = 0;
        psi_.psi_current = psi_.psi;
#endif
    }

    void store(int point)
    {
#if defined(__CUDA)
        if (!active()) { return; }
        assert(point >= 0 && point < points_ && psi_.get_current_k() == point);
        detail::kpoint_check(cudaMemcpy(host_ + size_t(point) * block_size_, psi_.psi,
                               block_size_ * sizeof(T), cudaMemcpyDeviceToHost));
#endif
    }

    void restore()
    {
#if defined(__CUDA)
        if (!active()) { return; }
        const int point = psi_.get_current_k();
        const int band = psi_.get_current_b();
        psi_.resize(points_, bands_, basis_);
        detail::kpoint_check(cudaMemcpy(psi_.psi, host_, size_t(points_) * block_size_ * sizeof(T),
                               cudaMemcpyHostToDevice));
        psi_.current_k = point;
        psi_.current_b = band;
        psi_.current_nbasis = psi_.ngk != nullptr ? psi_.ngk[point] : basis_;
        psi_.psi_bias = size_t(band) * basis_;
        psi_.psi_current = psi_.psi + psi_.psi_bias;
        detail::kpoint_check(cudaFreeHost(host_));
        host_ = nullptr;
#endif
    }

  private:
    Psi<T, Device>& psi_;
    const int points_;
    const int bands_;
    const int basis_;
    const size_t block_size_;
    T* host_ = nullptr;
};

} // namespace psi

#endif
