#ifndef HSOLVER_DAVIDSON_ACTION_H
#define HSOLVER_DAVIDSON_ACTION_H

#include "module_base/module_device/memory_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <type_traits>
#if defined(__CUDA)
#include <base/macros/macros.h>
#endif

namespace hsolver
{

/// Store Hamiltonian actions on Device or in pinned host memory. Host-backed
/// products use a bounded device panel; each GEMM retains its full contraction
/// dimension, so streaming changes storage without truncating the subspace.
template <typename T, typename Device>
class DavidsonAction
{
  public:
    DavidsonAction(int rows, int columns, bool host_backed)
        : rows_(rows), columns_(columns), panel_columns_(columns)
    {
        assert(rows_ > 0 && columns_ > 0);
        const size_t elements = size_t(rows_) * columns_;
#if defined(__CUDA)
        if (host_backed && std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            const auto status = cudaHostAlloc(reinterpret_cast<void**>(&host_),
                                              elements * sizeof(T), cudaHostAllocDefault);
            if (status == cudaErrorMemoryAllocation)
            {
                cudaGetLastError();
                host_ = nullptr;
            }
            else
            {
                cudaErrcheck(status);
            }
        }
#endif
        panel_columns_ = host_ ? std::min(64, columns_) : columns_;
        workspace_size_ = host_ ? std::max(size_t(rows_) * panel_columns_, size_t(columns_)) : elements;
        base_device::memory::resize_memory_op<T, Device>()(nullptr, device_, workspace_size_, "DAV::hphi");
        base_device::memory::set_memory_op<T, Device>()(nullptr, device_, 0, workspace_size_);
    }

    ~DavidsonAction()
    {
        base_device::memory::delete_memory_op<T, Device>()(nullptr, device_);
#if defined(__CUDA)
        if (host_) { cudaErrcheck(cudaFreeHost(host_)); }
#endif
    }

    DavidsonAction(const DavidsonAction&) = delete;
    DavidsonAction& operator=(const DavidsonAction&) = delete;

    bool host_backed() const { return host_ != nullptr; }
    int panel_columns() const { return panel_columns_; }
    T* resident() const
    {
        assert(!host_backed());
        return device_;
    }

    template <typename Action>
    void apply(const Action& action, T* basis, int basis_columns, int first, int count)
    {
        assert(first >= 0 && first + count <= columns_);
        if (!host_backed())
        {
            action(device_ + size_t(first) * rows_, basis, basis_columns, rows_, first, first + count - 1);
            return;
        }
        for (int offset = first; offset < first + count; offset += panel_columns_)
        {
            const int width = std::min(panel_columns_, first + count - offset);
            action(device_, basis, basis_columns, rows_, offset, offset + width - 1);
            store(offset, width, device_);
        }
    }

    void store(int first, int count, const T* source)
    {
        assert(first >= 0 && first + count <= columns_);
#if defined(__CUDA)
        if (host_backed())
        {
            cudaErrcheck(cudaMemcpy(host_ + size_t(first) * rows_, source,
                                   size_t(count) * rows_ * sizeof(T), cudaMemcpyDeviceToHost));
            return;
        }
#endif
        base_device::memory::synchronize_memory_op<T, Device, Device>()(
            nullptr, nullptr, device_ + size_t(first) * rows_, source, size_t(count) * rows_);
    }

    const T* panel(int first, int count)
    {
        assert(first >= 0 && first + count <= columns_ && count <= panel_columns_);
#if defined(__CUDA)
        if (host_backed())
        {
            cudaErrcheck(cudaMemcpy(device_, host_ + size_t(first) * rows_,
                                   size_t(count) * rows_ * sizeof(T), cudaMemcpyHostToDevice));
            return device_;
        }
#endif
        return device_ + size_t(first) * rows_;
    }

    /// Compute Hphi[:, :inner] * coefficients into a strided device matrix.
    void multiply(int inner, int output_columns, const T* coefficients, int coefficient_leading,
                  T* output, int output_leading)
    {
        assert(inner > 0 && inner <= columns_);
        const T one(1), zero(0);
        const int tile_rows = host_backed() ? std::min(size_t(rows_), workspace_size_ / inner) : rows_;
        for (int first = 0; first < rows_; first += tile_rows)
        {
            const int height = std::min(tile_rows, rows_ - first);
#if defined(__CUDA)
            if (host_backed())
            {
                cudaErrcheck(cudaMemcpy2D(device_, size_t(height) * sizeof(T),
                                         host_ + first, size_t(rows_) * sizeof(T),
                                         size_t(height) * sizeof(T), inner, cudaMemcpyHostToDevice));
            }
#endif
            gemm_op<T, Device>()(nullptr, 'N', 'N', height, output_columns, inner,
                                &one, device_, host_backed() ? height : rows_,
                                coefficients, coefficient_leading, &zero,
                                output + first, output_leading);
        }
    }

  private:
    const int rows_, columns_;
    int panel_columns_;
    size_t workspace_size_ = 0;
    T* device_ = nullptr;
    T* host_ = nullptr;
};

} // namespace hsolver

#endif
