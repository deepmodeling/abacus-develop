#ifndef DEVICE_CONTEXT_H_
#define DEVICE_CONTEXT_H_

#include <mutex>

#ifdef __MPI
#include "mpi.h"
#endif

namespace base_device {

/**
 * @brief Singleton class to manage GPU device context and initialization.
 *
 * This class provides a centralized way to:
 * 1. Initialize GPU device binding (only once)
 * 2. Query GPU device state (device_id, device_count, etc.)
 * 3. Ensure thread-safe initialization
 *
 * Usage:
 *   // Initialize (call once after MPI init and after determining device=gpu)
 *   DeviceContext::instance().init(MPI_COMM_WORLD);
 *
 *   // Query device info
 *   int dev_id = DeviceContext::instance().get_device_id();
 */
class DeviceContext {
public:
    /**
     * @brief Get the singleton instance of DeviceContext
     * @return Reference to the singleton instance
     */
    static DeviceContext& instance();

    /**
     * @brief Initialize GPU device binding.
     *
     * This function:
     * 1. Gets the local rank within the node using MPI_COMM_TYPE_SHARED
     * 2. Queries the number of available GPU devices
     * 3. Binds the current process to a GPU device (local_rank % device_count)
     *
     * @param comm MPI communicator (default: MPI_COMM_WORLD)
     * @note This function should only be called once. Subsequent calls are no-ops.
     * @note This function should only be called when device=gpu is confirmed.
     */
#ifdef __MPI
    void init(MPI_Comm comm = MPI_COMM_WORLD);
#else
    void init();
#endif

    /**
     * @brief Check if the DeviceContext has been initialized
     * @return true if init() has been called successfully
     */
    bool is_initialized() const { return initialized_; }

    /**
     * @brief Check if GPU is enabled and available
     * @return true if GPU device is bound and usable
     */
    bool is_gpu_enabled() const { return gpu_enabled_; }

    /**
     * @brief Get the bound GPU device ID
     * @return Device ID (0-based), or -1 if not initialized
     */
    int get_device_id() const { return device_id_; }

    /**
     * @brief Get the total number of GPU devices on this node
     * @return Number of GPU devices, or 0 if not initialized
     */
    int get_device_count() const { return device_count_; }

    /**
     * @brief Get the local MPI rank within the node
     * @return Local rank, or 0 if not initialized
     */
    int get_local_rank() const { return local_rank_; }

    // Disable copy and assignment
    DeviceContext(const DeviceContext&) = delete;
    DeviceContext& operator=(const DeviceContext&) = delete;

private:
    DeviceContext() = default;
    ~DeviceContext() = default;

    bool initialized_ = false;
    bool gpu_enabled_ = false;
    int device_id_ = -1;
    int device_count_ = 0;
    int local_rank_ = 0;

    std::mutex init_mutex_;
};

} // namespace base_device

#endif // DEVICE_CONTEXT_H_
