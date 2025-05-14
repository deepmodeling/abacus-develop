#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_TDDFT_EVOLVE_ELEC_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_TDDFT_EVOLVE_ELEC_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/module_container/ATen/core/tensor.h"     // ct::Tensor
#include "module_base/module_container/ATen/core/tensor_map.h" // TensorMap
#include "module_base/module_device/device.h"                  // base_device
#include "module_base/module_device/memory_op.h"               // memory operations
#include "module_base/scalapack_connector.h"                   // Cpxgemr2d
#include "module_esolver/esolver_ks_lcao.h"
#include "module_esolver/esolver_ks_lcao_tddft.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_psi/psi.h"

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
    std::cout << "Shape of " << name << ": [";
    for (int i = 0; i < tensor.shape().ndim(); ++i)
    {
        std::cout << tensor.shape().dim_size(i);
        if (i < tensor.shape().ndim() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

// Recursive print function
template <typename T>
inline void print_tensor_data_recursive(const T* data,
                                        const std::vector<int64_t>& shape,
                                        const std::vector<int64_t>& strides,
                                        int dim,
                                        std::vector<int64_t>& indices,
                                        const std::string& name)
{
    if (dim == shape.size())
    {
        // Recursion base case: print data when reaching the innermost dimension
        std::cout << name;
        for (size_t i = 0; i < indices.size(); ++i)
        {
            std::cout << "[" << indices[i] << "]";
        }
        std::cout << " = " << *data << std::endl;
        return;
    }
    // Recursively process the current dimension
    for (int64_t i = 0; i < shape[dim]; ++i)
    {
        indices[dim] = i;
        print_tensor_data_recursive(data + i * strides[dim], shape, strides, dim + 1, indices, name);
    }
}

// Generic print function
template <typename T>
inline void print_tensor_data(const ct::Tensor& tensor, const std::string& name)
{
    const std::vector<int64_t>& shape = tensor.shape().dims();
    const std::vector<int64_t>& strides = tensor.shape().strides();
    const T* data = tensor.data<T>();
    std::vector<int64_t> indices(shape.size(), 0);
    print_tensor_data_recursive(data, shape, strides, 0, indices, name);
}

// Specialization for std::complex<double>
template <>
inline void print_tensor_data<std::complex<double>>(const ct::Tensor& tensor, const std::string& name)
{
    const std::vector<int64_t>& shape = tensor.shape().dims();
    const std::vector<int64_t>& strides = tensor.shape().strides();
    const std::complex<double>* data = tensor.data<std::complex<double>>();
    std::vector<int64_t> indices(shape.size(), 0);
    print_tensor_data_recursive(data, shape, strides, 0, indices, name);
}

//------------------------ Debugging utility function ------------------------//

namespace module_tddft
{
#ifdef __MPI
//------------------------ MPI gathering and distributing functions ------------------------//
template <typename T>
void gatherPsi(const int myid,
               const int root_proc,
               T* psi_l,
               const Parallel_Orbitals& para_orb,
               ModuleESolver::Matrix_g<T>& psi_g)
{
    const int* desc_psi = para_orb.desc_wfc; // Obtain the descriptor from Parallel_Orbitals
    int ctxt = desc_psi[1];                  // BLACS context
    int nrows = desc_psi[2];                 // Global matrix row number
    int ncols = desc_psi[3];                 // Global matrix column number

    if (myid == root_proc)
    {
        psi_g.p.reset(new T[nrows * ncols]); // No need to delete[] since it is a shared_ptr
    }
    else
    {
        psi_g.p.reset(new T[nrows * ncols]); // Placeholder for non-root processes
    }

    // Set the descriptor of the global psi
    psi_g.desc.reset(new int[9]{1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows});
    psi_g.row = nrows;
    psi_g.col = ncols;

    // Call the Cpxgemr2d function in ScaLAPACK to collect the matrix data
    Cpxgemr2d(nrows, ncols, psi_l, 1, 1, const_cast<int*>(desc_psi), psi_g.p.get(), 1, 1, psi_g.desc.get(), ctxt);
}

template <typename T>
void distributePsi(const Parallel_Orbitals& para_orb, T* psi_l, const ModuleESolver::Matrix_g<T>& psi_g)
{
    const int* desc_psi = para_orb.desc_wfc; // Obtain the descriptor from Parallel_Orbitals
    int ctxt = desc_psi[1];                  // BLACS context
    int nrows = desc_psi[2];                 // Global matrix row number
    int ncols = desc_psi[3];                 // Global matrix column number

    // Call the Cpxgemr2d function in ScaLAPACK to distribute the matrix data
    Cpxgemr2d(nrows, ncols, psi_g.p.get(), 1, 1, psi_g.desc.get(), psi_l, 1, 1, const_cast<int*>(desc_psi), ctxt);
}
//------------------------ MPI gathering and distributing functions ------------------------//
#endif // __MPI

template <typename Device = base_device::DEVICE_CPU>
class Evolve_elec
{
    friend class ModuleESolver::ESolver_KS_LCAO<std::complex<double>, double>;

    // Template parameter is needed for the friend class declaration
    friend class ModuleESolver::ESolver_KS_LCAO_TDDFT<Device>;

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
                          std::complex<double>** Hk_laststep,
                          std::complex<double>** Sk_laststep,
                          ModuleBase::matrix& ekb,
                          std::ofstream& ofs_running,
                          const int htype,
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
} // namespace module_tddft
#endif
