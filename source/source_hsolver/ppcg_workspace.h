#ifndef HSOLVER_PPCG_WORKSPACE_H
#define HSOLVER_PPCG_WORKSPACE_H

#include "source_base/module_device/memory_op.h"
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hsolver
{

// Device-only storage deliberately has no host iterators or element access.
// Reallocation preserves existing columns using device-to-device copies.
template <typename T>
class PPCGDeviceBuffer
{
  public:
    PPCGDeviceBuffer() = default;
    PPCGDeviceBuffer(const size_t count, const T value) { assign(count, value); }
    ~PPCGDeviceBuffer() { base_device::memory::delete_memory_op<T, base_device::DEVICE_GPU>()(data_); }
    PPCGDeviceBuffer(const PPCGDeviceBuffer &) = delete;
    PPCGDeviceBuffer &operator=(const PPCGDeviceBuffer &) = delete;
    PPCGDeviceBuffer(PPCGDeviceBuffer &&other) noexcept { swap(other); }
    PPCGDeviceBuffer &operator=(PPCGDeviceBuffer &&other) noexcept
    {
        swap(other);
        return *this;
    }
    T *data() { return data_; }
    const T *data() const { return data_; }
    size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }
    void reserve(const size_t count)
    {
        if (count <= capacity_)
        {
            return;
        }
        T *replacement = nullptr;
        base_device::memory::resize_memory_op<T, base_device::DEVICE_GPU>()(replacement, count, "PPCG workspace");
        if (size_ > 0)
        {
            base_device::memory::synchronize_memory_op<T, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
                replacement, data_, size_);
        }
        base_device::memory::delete_memory_op<T, base_device::DEVICE_GPU>()(data_);
        data_ = replacement;
        capacity_ = count;
    }
    void resize(const size_t count)
    {
        reserve(count);
        if (count > size_)
        {
            T *tail = data_ + size_;
            const size_t added = count - size_;
            base_device::memory::set_memory_op<T, base_device::DEVICE_GPU>()(tail, 0, added);
        }
        size_ = count;
    }
    void resize(const size_t count, const T value)
    {
        require_zero(value);
        resize(count);
    }
    void assign(const size_t count, const T value)
    {
        require_zero(value);
        reserve(count);
        size_ = count;
        if (count > 0)
        {
            base_device::memory::set_memory_op<T, base_device::DEVICE_GPU>()(data_, 0, count);
        }
    }

  private:
    static void require_zero(const T value)
    {
        if (value != T(0))
        {
            throw std::invalid_argument("PPCG device workspaces support zero initialization only.");
        }
    }
    void swap(PPCGDeviceBuffer &other) noexcept
    {
        std::swap(data_, other.data_);
        std::swap(size_, other.size_);
        std::swap(capacity_, other.capacity_);
    }
    T *data_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
};

template <typename T, typename Device>
struct PPCGWorkspaceType
{
    using type = std::vector<T>;
};

template <typename T>
struct PPCGWorkspaceType<T, base_device::DEVICE_GPU>
{
    using type = PPCGDeviceBuffer<T>;
};

} // namespace hsolver
#endif
