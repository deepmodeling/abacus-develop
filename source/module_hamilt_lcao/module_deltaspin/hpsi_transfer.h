#ifndef DELTASPIN_HPSI_TRANSFER_H
#define DELTASPIN_HPSI_TRANSFER_H

#include <cassert>
#include <complex>
#include <base/macros/macros.h>

namespace spinconstrain
{

/// Restore immutable pinned-host H psi blocks on the default CUDA stream.
/// One host thread consumes the blocks; the host storage outlives this object.
/// Direct uploads retain no duplicate device wavefunction block. Destruction
/// waits for pending uploads before the caller can release the host storage.
class HPsiTransfer
{
    using Complex = std::complex<double>;
    const Complex* host_;
    const int points_, bands_, source_leading_;
    cudaEvent_t uploaded_ = nullptr;

  public:
    HPsiTransfer(const Complex* host, int points, int bands, int source_leading)
        : host_(host), points_(points), bands_(bands), source_leading_(source_leading)
    {
        assert(host != nullptr && points > 0 && bands > 0 && source_leading > 0);
    }

    HPsiTransfer(const HPsiTransfer&) = delete;
    HPsiTransfer& operator=(const HPsiTransfer&) = delete;

    ~HPsiTransfer()
    {
        if (uploaded_ == nullptr) { return; }
        cudaErrcheck(cudaEventSynchronize(uploaded_));
        cudaErrcheck(cudaEventDestroy(uploaded_));
    }

    /// Copy bands columns with source_leading rows into compact output columns.
    /// Requests may repeat or reorder k points and vary the output row count.
    void restore(int ik, Complex* output, int leading)
    {
        assert(ik >= 0 && ik < points_ && output != nullptr);
        assert(leading > 0 && leading <= source_leading_);
        if (uploaded_ == nullptr)
        {
            cudaErrcheck(cudaEventCreateWithFlags(&uploaded_, cudaEventDisableTiming));
        }
        const size_t pitch = size_t(leading) * sizeof(Complex);
        const size_t source_pitch = size_t(source_leading_) * sizeof(Complex);
        const size_t offset = size_t(ik) * bands_ * source_leading_;
        cudaErrcheck(cudaMemcpy2DAsync(output, pitch, host_ + offset, source_pitch,
            pitch, bands_, cudaMemcpyHostToDevice, nullptr));
        // All uploads share the consumer stream. Its last recorded event
        // therefore covers every host read, including repeated requests.
        cudaErrcheck(cudaEventRecord(uploaded_, nullptr));
    }
};

} // namespace spinconstrain

#endif
