/**
 * @file deltaspin_pw_cache.h
 * @brief PW-basis subspace data cache for DeltaSpin, decoupled from SpinConstrain.
 *
 * @par Purpose
 * In the PW basis, the subspace Hamiltonian H_sub = <psi|H|psi>, overlap S_sub
 * and becp coefficients are expensive to compute. They are cached on the first
 * cal_mw_from_lambda() call and reused across multiple lambda steps within the
 * same SCF iteration, then freed after the final subspace diagonalization in
 * update_psi_charge_pw_{cpu,gpu}().
 *
 * This class owns the three host buffers and raw device pointers plus the lambda
 * snapshot taken when the cache was filled. It encapsulates the CPU vs GPU
 * allocation/free difference behind allocate()/release(), while still exposing
 * raw per-k pointers (h_k/s_k/becp_k) because the hsolver subspace routines and
 * GPU memcpy ops require raw pointers.
 *
 * @par Layout (same as before, unchanged)
 * - h(ik)[i * nbands + j]: H_sub for k-point ik
 * - s(ik): same layout for overlap S_sub
 * - becp(ik)[ib * nkb * npol + ip]: becp coefficients
 */
#ifndef DELTASPIN_PW_CACHE_H
#define DELTASPIN_PW_CACHE_H

#include <complex>
#include <vector>

#include "source_base/vector3.h"
#include "source_base/module_device/memory_op.h"

namespace spinconstrain
{
namespace pw
{

/**
 * @brief Owning cache of PW subspace H/S/becp data plus the lambda snapshot.
 *
 * The buffer element type is std::complex<double> because the PW DeltaSpin path
 * is always instantiated on complex wavefunctions; the legacy TK=double stub
 * never allocates it.
 */
class SubspaceCache
{
  public:
    SubspaceCache() = default;

    // Owns cache memory; non-copyable, non-movable to keep ownership unambiguous.
    SubspaceCache(const SubspaceCache&) = delete;
    SubspaceCache& operator=(const SubspaceCache&) = delete;

    /// True when the subspace buffers are allocated.
    bool allocated() const
    {
#if ((defined __CUDA) || (defined __ROCM))
        return !sub_h_save_.empty() || sub_h_save_gpu_ != nullptr;
#else
        return !sub_h_save_.empty();
#endif
    }

    /// Lambda values captured when the cache was filled.
    std::vector<ModuleBase::Vector3<double>>& lambda_in_sub() { return lambda_in_sub_; }
    const std::vector<ModuleBase::Vector3<double>>& lambda_in_sub() const { return lambda_in_sub_; }

    /// Raw base pointers (needed by hsolver subspace ops and GPU memcpy).
    std::complex<double>* h()
    {
#if ((defined __CUDA) || (defined __ROCM))
        if (sub_h_save_gpu_ != nullptr)
        {
            return sub_h_save_gpu_;
        }
#endif
        return sub_h_save_.data();
    }
    std::complex<double>* s()
    {
#if ((defined __CUDA) || (defined __ROCM))
        if (sub_s_save_gpu_ != nullptr)
        {
            return sub_s_save_gpu_;
        }
#endif
        return sub_s_save_.data();
    }
    std::complex<double>* becp()
    {
#if ((defined __CUDA) || (defined __ROCM))
        if (becp_save_gpu_ != nullptr)
        {
            return becp_save_gpu_;
        }
#endif
        return becp_save_.data();
    }

    /// Per-k-point views.
    std::complex<double>* h_k(int ik, int nbands)
    {
        std::complex<double>* base = h();
        return base == nullptr ? nullptr : base + ik * nbands * nbands;
    }
    std::complex<double>* s_k(int ik, int nbands)
    {
        std::complex<double>* base = s();
        return base == nullptr ? nullptr : base + ik * nbands * nbands;
    }
    std::complex<double>* becp_k(int ik, int size_becp)
    {
        std::complex<double>* base = becp();
        return base == nullptr ? nullptr : base + ik * size_becp;
    }

    /**
     * @brief Allocate the three buffers on the host (CPU path) with resize().
     * No-op if already allocated.
     */
    void allocate_cpu(int nbands, int nk, int size_becp)
    {
        if (allocated())
        {
            return;
        }
        sub_h_save_.resize(nbands * nbands * nk);
        sub_s_save_.resize(nbands * nbands * nk);
        becp_save_.resize(size_becp * nk);
    }

    /**
     * @brief Release the host (CPU) buffers and their storage capacity.
     */
    void release_cpu()
    {
        std::vector<std::complex<double>>().swap(sub_h_save_);
        std::vector<std::complex<double>>().swap(sub_s_save_);
        std::vector<std::complex<double>>().swap(becp_save_);
    }

#if ((defined __CUDA) || (defined __ROCM))
    /**
     * @brief Allocate the three buffers on the device (GPU path).
     * No-op if already allocated.
     */
    void allocate_gpu(int nbands, int nk, int size_becp)
    {
        if (allocated())
        {
            return;
        }
        using mem = base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
        mem()(sub_h_save_gpu_, nbands * nbands * nk);
        mem()(sub_s_save_gpu_, nbands * nbands * nk);
        mem()(becp_save_gpu_, size_becp * nk);
    }

    /**
     * @brief Release the device (GPU) buffers.
     */
    void release_gpu()
    {
        using del = base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
        del()(sub_h_save_gpu_);
        del()(sub_s_save_gpu_);
        del()(becp_save_gpu_);
        sub_h_save_gpu_ = nullptr;
        sub_s_save_gpu_ = nullptr;
        becp_save_gpu_ = nullptr;
    }
#endif // __CUDA || __ROCM

  private:
    std::vector<std::complex<double>> sub_h_save_; ///< Host subspace Hamiltonian for all k-points
    std::vector<std::complex<double>> sub_s_save_; ///< Host subspace overlap matrix for all k-points
    std::vector<std::complex<double>> becp_save_;  ///< Host becp coefficients for all k-points
#if ((defined __CUDA) || (defined __ROCM))
    std::complex<double>* sub_h_save_gpu_ = nullptr; ///< Device subspace Hamiltonian for all k-points
    std::complex<double>* sub_s_save_gpu_ = nullptr; ///< Device subspace overlap matrix for all k-points
    std::complex<double>* becp_save_gpu_ = nullptr;  ///< Device becp coefficients for all k-points
#endif
    std::vector<ModuleBase::Vector3<double>> lambda_in_sub_; ///< Lambda when the cache was saved
};

} // namespace pw
} // namespace spinconstrain

#endif // DELTASPIN_PW_CACHE_H
