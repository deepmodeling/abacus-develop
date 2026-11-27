#ifndef HSOLVER_LOBPCG_INTERNAL_H_
#define HSOLVER_LOBPCG_INTERNAL_H_

// Shared implementation details, included only by the LOBPCG .cpp files.
#include "diago_lobpcg.h"
#include "source_base/module_container/ATen/core/tensor_map.h"
#include "source_base/parallel_device.h"

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace hsolver
{
namespace lobpcg_internal
{
/// Allocate and fill a host Tensor for small projected matrices or scalar work arrays.
template <typename T>
ct::Tensor host_tensor(const int64_t count, const T value)
{
    const ct::DataType type = ct::DataTypeToEnum<T>::value;
    ct::Tensor tensor(type, ct::DeviceType::CpuDevice, {count});
    std::fill_n(tensor.data<T>(), count, value);
    return tensor;
}

/// Borrow host storage when types match; otherwise own a converted copy. The view must not outlive input.
template <typename T, typename U>
ct::Tensor host_tensor_cast_or_view(U* input, const int64_t count)
{
    const ct::DataType type = ct::DataTypeToEnum<T>::value;
    if (std::is_same<T, U>::value)
    {
        ct::TensorMap view(input, type, ct::DeviceType::CpuDevice, {count});
        return ct::Tensor(std::move(view));
    }
    ct::Tensor tensor(type, ct::DeviceType::CpuDevice, {count});
    std::copy_n(input, count, tensor.data<T>());
    return tensor;
}

} // namespace lobpcg_internal

/// Share root-owned coefficients or status within the eigensolve communicator.
template <typename T, typename Device>
template <typename Value>
void DiagoLOBPCG<T, Device>::broadcast(Value* data, const int count) const
{
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::bcast_data(data, count, comm_.comm, 0);
    }
#endif
}

} // namespace hsolver

#endif // HSOLVER_LOBPCG_INTERNAL_H_
