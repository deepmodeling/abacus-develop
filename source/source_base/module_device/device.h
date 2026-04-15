#ifndef MODULE_DEVICE_H_
#define MODULE_DEVICE_H_

#include "types.h"
#include "device_helpers.h"
#include <fstream>
#include <mutex>

#ifdef __MPI
#include "mpi.h"
#endif

namespace base_device
{

namespace information
{

/**
 * @brief Get the device name
 * for source_esolver
 */
std::string get_device_name(std::string device_flag);

/**
 * @brief Get the device number
 * for source_esolver
 */
int get_device_num(std::string device_flag);

/**
 * @brief Output the device information
 * for source_esolver
 */
void output_device_info(std::ostream& output);

/**
 * @brief Safely probes for GPU availability without exiting on error.
 * @return True if at least one GPU is found and usable, false otherwise.
 */
bool probe_gpu_availability();

/**
 * @brief Get the device flag object
 * for source_io PARAM.inp.device
 */
std::string get_device_flag(const std::string& device,
                            const std::string& basis_type);

#if __MPI
/**
 * @brief Get the local rank within the node using MPI_COMM_TYPE_SHARED
 * @param mpi_comm MPI communicator (default: MPI_COMM_WORLD)
 * @return Local rank within the node
 */
int get_node_rank_with_mpi_shared(const MPI_Comm mpi_comm = MPI_COMM_WORLD);
#endif

template <typename Device>
void print_device_info(const Device* dev, std::ofstream& ofs_device)
{
    return;
}

template <typename Device>
void record_device_memory(const Device* dev, std::ofstream& ofs_device, std::string str, size_t size)
{
    return;
}

#if defined(__CUDA) || defined(__ROCM)
template <>
void print_device_info<base_device::DEVICE_GPU>(const base_device::DEVICE_GPU *ctx, std::ofstream &ofs_device);

template <>
void record_device_memory<base_device::DEVICE_GPU>(const base_device::DEVICE_GPU* dev, std::ofstream& ofs_device, std::string str, size_t size);
#endif

} // end of namespace information

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
     * 1. Gets the local rank within the node using MPI_COMM_TYPE_SHARED (MPI_COMM_WORLD)
     * 2. Queries the number of available GPU devices
     * 3. Binds the current process to a GPU device (local_rank % device_count)
     *
     * @note This function should only be called once. Subsequent calls are no-ops.
     * @note This function should only be called when device=gpu is confirmed.
     * @note In MPI builds, uses MPI_COMM_WORLD internally.
     */
    void init();

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

    /**
     * @brief Set the device type (CpuDevice, GpuDevice, or DspDevice)
     * @param type The device type
     */
    void set_device_type(AbacusDevice_t type) { device_type_ = type; }

    /**
     * @brief Get the device type
     * @return AbacusDevice_t The device type
     */
    AbacusDevice_t get_device_type() const { return device_type_; }

    /**
     * @brief Check if the device is CPU
     * @return true if the device is CPU
     */
    bool is_cpu() const { return device_type_ == CpuDevice; }

    /**
     * @brief Check if the device is GPU
     * @return true if the device is GPU
     */
    bool is_gpu() const { return device_type_ == GpuDevice; }

    /**
     * @brief Check if the device is DSP
     * @return true if the device is DSP
     */
    bool is_dsp() const { return device_type_ == DspDevice; }

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
    AbacusDevice_t device_type_ = CpuDevice;

    std::mutex init_mutex_;
};

/**
 * @brief Get the device type enum from DeviceContext (runtime version).
 * @param ctx Pointer to DeviceContext
 * @return AbacusDevice_t enum value
 */
inline AbacusDevice_t get_device_type(const DeviceContext* ctx)
{
    return ctx->get_device_type();
}

} // end of namespace base_device

#endif // MODULE_DEVICE_H_
