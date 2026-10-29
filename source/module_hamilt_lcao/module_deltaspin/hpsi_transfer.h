#ifndef DELTASPIN_HPSI_TRANSFER_H
#define DELTASPIN_HPSI_TRANSFER_H

#include <cassert>
#include <complex>
#include <base/macros/macros.h>

namespace spinconstrain
{

/// Restore immutable pinned-host H psi blocks on the default CUDA stream.
/// One host thread consumes the blocks; the host storage outlives this object.
/// A single device block prefetches the next k point on a nonblocking stream.
class HPsiTransfer
{
    using Complex = std::complex<double>;
    const Complex* host_;
    const int points_, bands_, source_leading_;
    Complex* storage_ = nullptr;
    cudaStream_t stream_ = nullptr;
    cudaEvent_t ready_ = nullptr, consumed_ = nullptr;
    int cached_index_ = -1;
    bool has_consumer_ = false;
    bool allocation_attempted_;

    void allocate()
    {
        // Davidson arrays and the current projectors take priority over this
        // optional buffer. A failed allocation selects direct uploads for the
        // remainder of this solve, without repeating the allocation attempt.
        allocation_attempted_ = true;
        const auto status = cudaMalloc(reinterpret_cast<void**>(&storage_),
            size_t(bands_) * source_leading_ * sizeof(Complex));
        if (status == cudaErrorMemoryAllocation)
        {
            storage_ = nullptr;
            cudaGetLastError();
            return;
        }
        cudaErrcheck(status);
        cudaErrcheck(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
        cudaErrcheck(cudaEventCreateWithFlags(&ready_, cudaEventDisableTiming));
        cudaErrcheck(cudaEventCreateWithFlags(&consumed_, cudaEventDisableTiming));
    }

    void prepare(int ik)
    {
        // The copy stream cannot overwrite the block until the preceding
        // device-to-device restore has consumed it on the default stream.
        if (has_consumer_) { cudaErrcheck(cudaStreamWaitEvent(stream_, consumed_, 0)); }
        const size_t count = size_t(bands_) * source_leading_;
        cudaErrcheck(cudaMemcpyAsync(storage_, host_ + size_t(ik) * count,
            count * sizeof(Complex), cudaMemcpyHostToDevice, stream_));
        cudaErrcheck(cudaEventRecord(ready_, stream_));
        cached_index_ = ik;
    }

  public:
    HPsiTransfer(const Complex* host, int points, int bands, int source_leading)
        : host_(host), points_(points), bands_(bands), source_leading_(source_leading),
          allocation_attempted_(points == 1)
    {
        assert(host != nullptr && points > 0 && bands > 0 && source_leading > 0);
    }

    HPsiTransfer(const HPsiTransfer&) = delete;
    HPsiTransfer& operator=(const HPsiTransfer&) = delete;

    ~HPsiTransfer()
    {
        if (storage_ == nullptr) { return; }
        cudaErrcheck(cudaStreamSynchronize(stream_));
        if (has_consumer_) { cudaErrcheck(cudaEventSynchronize(consumed_)); }
        cudaErrcheck(cudaFree(storage_));
        cudaErrcheck(cudaEventDestroy(ready_));
        cudaErrcheck(cudaEventDestroy(consumed_));
        cudaErrcheck(cudaStreamDestroy(stream_));
    }

    /// Copy bands columns with source_leading rows into compact output columns.
    /// Reordered or repeated k points reload the requested block when necessary.
    void restore(int ik, Complex* output, int leading)
    {
        assert(ik >= 0 && ik < points_ && output != nullptr);
        assert(leading > 0 && leading <= source_leading_);
        if (!allocation_attempted_) { allocate(); }
        const size_t pitch = size_t(leading) * sizeof(Complex);
        const size_t source_pitch = size_t(source_leading_) * sizeof(Complex);
        if (storage_ == nullptr)
        {
            const size_t offset = size_t(ik) * bands_ * source_leading_;
            cudaErrcheck(cudaMemcpy2D(output, pitch, host_ + offset, source_pitch,
                pitch, bands_, cudaMemcpyHostToDevice));
            return;
        }
        if (cached_index_ != ik) { prepare(ik); }
        cudaErrcheck(cudaStreamWaitEvent(nullptr, ready_, 0));
        cudaErrcheck(cudaMemcpy2DAsync(output, pitch, storage_, source_pitch,
            pitch, bands_, cudaMemcpyDeviceToDevice, nullptr));
        cudaErrcheck(cudaEventRecord(consumed_, nullptr));
        has_consumer_ = true;
        if (ik + 1 < points_) { prepare(ik + 1); }
    }
};

} // namespace spinconstrain

#endif
