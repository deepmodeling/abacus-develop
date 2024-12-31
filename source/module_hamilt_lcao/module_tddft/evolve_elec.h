#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_TDDFT_EVOLVE_ELEC_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_TDDFT_EVOLVE_ELEC_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/module_container/ATen/core/tensor.h"     // container::Tensor
#include "module_base/module_container/ATen/core/tensor_map.h" // TensorMap
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
inline void print_tensor_shape(const container::Tensor& tensor, const std::string& name)
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
inline void print_tensor_data(const container::Tensor& tensor, const std::string& name)
{
    const std::vector<int64_t>& shape = tensor.shape().dims();
    const std::vector<int64_t>& strides = tensor.shape().strides();
    const T* data = tensor.data<T>();
    std::vector<int64_t> indices(shape.size(), 0);
    print_tensor_data_recursive(data, shape, strides, 0, indices, name);
}

// Specialization for std::complex<double>
template <>
inline void print_tensor_data<std::complex<double>>(const container::Tensor& tensor, const std::string& name)
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
class Evolve_elec
{

    friend class ELEC_scf;
    friend class ModuleESolver::ESolver_KS_LCAO<std::complex<double>, double>;
    friend class ModuleESolver::ESolver_KS_LCAO_TDDFT;

  public:
    Evolve_elec();
    ~Evolve_elec();

    static double td_force_dt;
    static bool td_vext;
    static std::vector<int> td_vext_dire_case;
    // Output dipole, efield, current or not
    static bool out_dipole;
    static bool out_efield;

    static double td_print_eij; // the threshold to output Eij elements
    static int td_edm;          // 0: new edm method   1: old edm method

  private:
    static void solve_psi(const int& istep,
                          const int nband,
                          const int nlocal,
                          hamilt::Hamilt<std::complex<double>>* phm,
                          Parallel_Orbitals& para_orb,
                          psi::Psi<std::complex<double>>* psi,
                          psi::Psi<std::complex<double>>* psi_laststep,
                          std::complex<double>** Hk_laststep,
                          std::complex<double>** Sk_laststep,
                          ModuleBase::matrix& ekb,
                          int htype,
                          int propagator,
                          const int& nks);
};
} // namespace module_tddft
#endif
