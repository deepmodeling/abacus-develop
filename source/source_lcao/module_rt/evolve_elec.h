#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_TDDFT_EVOLVE_ELEC_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_TDDFT_EVOLVE_ELEC_H

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/module_container/ATen/core/tensor.h"     // ct::Tensor
#include "source_base/module_container/ATen/core/tensor_map.h" // TensorMap
#include "source_base/module_device/device.h"                  // base_device
#include "source_base/module_device/memory_op.h"               // memory operations
#include "source_base/module_external/scalapack_connector.h"   // Cpxgemr2d
#include "source_esolver/esolver_ks_lcao.h"
#include "source_esolver/esolver_ks_lcao_tddft.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_rt/gather_mat.h" // MPI gathering and distributing functions
#include "source_lcao/module_rt/kernels/cublasmp_context.h"
#include "source_psi/psi.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to evolve the electronic wave functions
// in TDDFT in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

//------------------------ Debugging utility function ------------------------//

// Print the shape of a Tensor
inline void print_tensor_shape(const ct::Tensor& tensor, const std::string& name)
{
    GlobalV::ofs_running << "Shape of " << name << ": [";
    for (int i = 0; i < tensor.shape().ndim(); ++i)
    {
        GlobalV::ofs_running << tensor.shape().dim_size(i);
        if (i < tensor.shape().ndim() - 1)
        {
            GlobalV::ofs_running << ", ";
        }
    }
    GlobalV::ofs_running << "]" << std::endl;
}

// Recursive print function
template <typename T>
inline void print_single_element(const T& val, double threshold)
{
    double clean_val = (std::abs(val) < threshold) ? 0.0 : static_cast<double>(val);
    GlobalV::ofs_running << std::fixed << std::setprecision(6) << clean_val;
}
inline void print_single_element(const std::complex<double>& val, double threshold)
{
    double re = (std::abs(val.real()) < threshold) ? 0.0 : val.real();
    double im = (std::abs(val.imag()) < threshold) ? 0.0 : val.imag();
    GlobalV::ofs_running << std::fixed << std::setprecision(6) << "(" << re << "," << im << ")";
}

template <typename T>
inline void print_tensor_data_recursive(const T* data,
                                        const std::vector<int64_t>& shape,
                                        const std::vector<int64_t>& strides,
                                        int dim,
                                        std::vector<int64_t>& indices,
                                        const std::string& name,
                                        const double threshold = 1e-10)
{
    if (dim == shape.size())
    {
        GlobalV::ofs_running << name;
        for (size_t i = 0; i < indices.size(); ++i)
        {
            GlobalV::ofs_running << "[" << indices[i] << "]";
        }
        GlobalV::ofs_running << " = ";

        print_single_element(*data, threshold);

        GlobalV::ofs_running << std::endl;
        return;
    }

    for (int64_t i = 0; i < shape[dim]; ++i)
    {
        indices[dim] = i;
        print_tensor_data_recursive(data + i * strides[dim], shape, strides, dim + 1, indices, name, threshold);
    }
}

template <typename T>
inline void print_tensor_data(const ct::Tensor& tensor, const std::string& name)
{
    const ct::Tensor* p_tensor = &tensor;
    ct::Tensor cpu_tensor_buffer;

    if (tensor.device_type() != ct::DeviceType::CpuDevice)
    {
        cpu_tensor_buffer = tensor.to_device<ct::DEVICE_CPU>();
        p_tensor = &cpu_tensor_buffer;
    }

    const std::vector<int64_t>& shape = p_tensor->shape().dims();
    const std::vector<int64_t>& strides = p_tensor->shape().strides();

    const T* data = p_tensor->data<T>();

    std::vector<int64_t> indices(shape.size(), 0);
    print_tensor_data_recursive(data, shape, strides, 0, indices, name);
}

template <>
inline void print_tensor_data<std::complex<double>>(const ct::Tensor& tensor, const std::string& name)
{
    const ct::Tensor* p_tensor = &tensor;
    ct::Tensor cpu_tensor_buffer;

    if (tensor.device_type() != ct::DeviceType::CpuDevice)
    {
        cpu_tensor_buffer = tensor.to_device<ct::DEVICE_CPU>();
        p_tensor = &cpu_tensor_buffer;
    }

    const std::vector<int64_t>& shape = p_tensor->shape().dims();
    const std::vector<int64_t>& strides = p_tensor->shape().strides();

    const std::complex<double>* data = p_tensor->data<std::complex<double>>();

    std::vector<int64_t> indices(shape.size(), 0);
    print_tensor_data_recursive(data, shape, strides, 0, indices, name);
}

//------------------------ Debugging utility function ------------------------//

namespace module_rt
{
template <typename Device = base_device::DEVICE_CPU>
class Evolve_elec
{
    friend class ModuleESolver::ESolver_KS_LCAO<std::complex<double>, double>;

    // Template parameter is needed for the friend class declaration
    friend class ModuleESolver::ESolver_KS_LCAO_TDDFT<double, Device>;
    friend class ModuleESolver::ESolver_KS_LCAO_TDDFT<std::complex<double>, Device>;

  public:
    Evolve_elec();
    ~Evolve_elec();

  private:
    static void solve_psi(const int& istep,
                          const int nband,
                          const int nlocal,
                          const int& nks,
                          hamilt::Hamilt<std::complex<double>>* phm,
                          Parallel_Orbitals& para_orb,
                          psi::Psi<std::complex<double>>* psi,
                          psi::Psi<std::complex<double>>* psi_laststep,
                          ct::Tensor& Hk_laststep,
                          ct::Tensor& Sk_laststep,
                          ModuleBase::matrix& ekb,
                          std::ofstream& ofs_running,
                          const int propagator,
                          const bool use_tensor,
                          const bool use_lapack);

    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    static ct::DeviceType ct_device_type;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Memory operations
    using syncmem_double_h2d_op = base_device::memory::synchronize_memory_op<double, Device, base_device::DEVICE_CPU>;
    using syncmem_double_d2h_op = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, Device>;
    using syncmem_complex_h2d_op
        = base_device::memory::synchronize_memory_op<std::complex<double>, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op
        = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, Device>;
};
} // namespace module_rt
#endif
