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
 * This class owns the three raw device/host pointers plus the lambda snapshot
 * taken when the cache was filled, replacing the ad-hoc new[]/delete[] that used
 * to live as public SpinConstrain members. It encapsulates the CPU vs GPU
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

    // Owns raw memory; non-copyable, non-movable to keep ownership unambiguous.
    SubspaceCache(const SubspaceCache&) = delete;
    SubspaceCache& operator=(const SubspaceCache&) = delete;

    /// True when the subspace buffers are allocated.
    bool allocated() const { return sub_h_save_ != nullptr; }

    /// Lambda values captured when the cache was filled.
    std::vector<ModuleBase::Vector3<double>>& lambda_in_sub() { return lambda_in_sub_; }
    const std::vector<ModuleBase::Vector3<double>>& lambda_in_sub() const { return lambda_in_sub_; }

    /// Raw base pointers (needed by hsolver subspace ops and GPU memcpy).
    std::complex<double>* h() { return sub_h_save_; }
    std::complex<double>* s() { return sub_s_save_; }
    std::complex<double>* becp() { return becp_save_; }

    /// Per-k-point views.
    std::complex<double>* h_k(int ik, int nbands) { return sub_h_save_ + ik * nbands * nbands; }
    std::complex<double>* s_k(int ik, int nbands) { return sub_s_save_ + ik * nbands * nbands; }
    std::complex<double>* becp_k(int ik, int size_becp) { return becp_save_ + ik * size_becp; }

    /**
     * @brief Allocate the three buffers on the host (CPU path) with new[].
     * No-op if already allocated.
     */
    void allocate_cpu(int nbands, int nk, int size_becp)
    {
        if (allocated())
        {
            return;
        }
        sub_h_save_ = new std::complex<double>[nbands * nbands * nk];
        sub_s_save_ = new std::complex<double>[nbands * nbands * nk];
        becp_save_ = new std::complex<double>[size_becp * nk];
    }

    /**
     * @brief Release the host (CPU) buffers with delete[].
     */
    void release_cpu()
    {
        delete[] sub_h_save_;
        delete[] sub_s_save_;
        delete[] becp_save_;
        sub_h_save_ = nullptr;
        sub_s_save_ = nullptr;
        becp_save_ = nullptr;
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
        mem()(sub_h_save_, nbands * nbands * nk);
        mem()(sub_s_save_, nbands * nbands * nk);
        mem()(becp_save_, size_becp * nk);
    }

    /**
     * @brief Release the device (GPU) buffers.
     */
    void release_gpu()
    {
        using del = base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>;
        del()(sub_h_save_);
        del()(sub_s_save_);
        del()(becp_save_);
        sub_h_save_ = nullptr;
        sub_s_save_ = nullptr;
        becp_save_ = nullptr;
    }
#endif // __CUDA || __ROCM

  private:
    std::complex<double>* sub_h_save_ = nullptr; ///< Cached subspace Hamiltonian for all k-points
    std::complex<double>* sub_s_save_ = nullptr; ///< Cached subspace overlap matrix for all k-points
    std::complex<double>* becp_save_ = nullptr;  ///< Cached becp coefficients for all k-points
    std::vector<ModuleBase::Vector3<double>> lambda_in_sub_; ///< Lambda when the cache was saved
};

} // namespace pw
} // namespace spinconstrain

#endif // DELTASPIN_PW_CACHE_H
