/**
 * @file propagator.h
 * @brief compute propagtor to evolve the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "module_base/constants.h"
#include "module_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "module_basis/module_ao/parallel_orbitals.h"

#include <complex>

namespace module_tddft
{
//--------------------------------- Utility function ---------------------------------//
#ifdef __MPI
inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}
#endif // __MPI

// Auxiliary function: process non-complex types, return value 1.0
template <typename T>
inline T init_value(typename std::enable_if<!std::is_same<T, std::complex<float>>::value
                                            && !std::is_same<T, std::complex<double>>::value>::type* = nullptr)
{
    return T(1.0);
}

// Auxiliary function: process complex types, return value 1.0 + 0.0i
template <typename T>
inline T init_value(typename std::enable_if<std::is_same<T, std::complex<float>>::value
                                            || std::is_same<T, std::complex<double>>::value>::type* = nullptr)
{
    return T(1.0, 0.0);
}

// Create an identity matrix of size n×n
template <typename T>
ct::Tensor create_identity_matrix(const int n, ct::DeviceType device = ct::DeviceType::CpuDevice)
{
    // Choose the data type of the Tensor
    ct::DataType data_type;
    if (std::is_same<T, float>::value)
    {
        data_type = ct::DataType::DT_FLOAT;
    }
    else if (std::is_same<T, double>::value)
    {
        data_type = ct::DataType::DT_DOUBLE;
    }
    else if (std::is_same<T, std::complex<float>>::value)
    {
        data_type = ct::DataType::DT_COMPLEX;
    }
    else if (std::is_same<T, std::complex<double>>::value)
    {
        data_type = ct::DataType::DT_COMPLEX_DOUBLE;
    }
    else
    {
        static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value
                          || std::is_same<T, std::complex<float>>::value
                          || std::is_same<T, std::complex<double>>::value,
                      "Unsupported data type!");
    }

    ct::Tensor tensor(data_type, device, ct::TensorShape({n, n}));
    tensor.zero();

    // Set the diagonal elements to 1
    if (device == ct::DeviceType::CpuDevice)
    {
        // For CPU, we can directly access the data
        T* data_ptr = tensor.data<T>();
        for (int i = 0; i < n; ++i)
        {
            data_ptr[i * n + i] = init_value<T>();
        }
    }
    else if (device == ct::DeviceType::GpuDevice)
    {
        // For GPU, we need to use a kernel to set the diagonal elements
        T* data_ptr = tensor.data<T>();
        for (int i = 0; i < n; ++i)
        {
            T value = init_value<T>();
            ct::kernels::set_memory<T, ct::DEVICE_GPU>()(data_ptr + i * n + i, value, 1);
        }
    }

    return tensor;
}
//--------------------------------- Utility function ---------------------------------//

class Propagator
{
  public:
    Propagator(const int ptype, const Parallel_Orbitals* pv, const double& dt)
    {
        this->ptype = ptype;
        this->ParaV = pv;
        this->dt = dt / ModuleBase::AU_to_FS;
    }
    ~Propagator();

#ifdef __MPI
    /**
     *  @brief compute propagator
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] H_laststep H(t)
     * @param[in] print_matirx print internal matrix or not
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator(const int nlocal,
                            const std::complex<double>* Stmp,
                            const std::complex<double>* Htmp,
                            const std::complex<double>* H_laststep,
                            std::complex<double>* U_operator,
                            std::ofstream& ofs_running,
                            const int print_matrix) const;

    template <typename Device>
    void compute_propagator_tensor(const int nlocal,
                                   const ct::Tensor& Stmp,
                                   const ct::Tensor& Htmp,
                                   const ct::Tensor& H_laststep,
                                   ct::Tensor& U_operator,
                                   std::ofstream& ofs_running,
                                   const int print_matrix,
                                   const bool use_lapack) const;
#endif // __MPI

  private:
    int ptype; // type of propagator
    const Parallel_Orbitals* ParaV;
    double dt; // time step

#ifdef __MPI

    /**
     *  @brief compute propagator of method Crank-Nicolson
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] print_matirx print internal matrix or not
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator_cn2(const int nlocal,
                                const std::complex<double>* Stmp,
                                const std::complex<double>* Htmp,
                                std::complex<double>* U_operator,
                                std::ofstream& ofs_running,
                                const int print_matrix) const;

    void compute_propagator_cn2_tensor(const int nlocal,
                                       const ct::Tensor& Stmp,
                                       const ct::Tensor& Htmp,
                                       ct::Tensor& U_operator,
                                       std::ofstream& ofs_running,
                                       const int print_matrix) const;

    template <typename Device>
    void compute_propagator_cn2_tensor_lapack(const int nlocal,
                                              const ct::Tensor& Stmp,
                                              const ct::Tensor& Htmp,
                                              ct::Tensor& U_operator,
                                              std::ofstream& ofs_running,
                                              const int print_matrix) const;

    /**
     *  @brief compute propagator of method 4th Taylor
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] print_matirx print internal matrix or not
     * @param[in] tag a parametre different for 4th Taylor and ETRS
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator_taylor(const int nlocal,
                                   const std::complex<double>* Stmp,
                                   const std::complex<double>* Htmp,
                                   std::complex<double>* U_operator,
                                   std::ofstream& ofs_running,
                                   const int print_matrix,
                                   const int tag) const;

    /**
     *  @brief compute propagator of method ETRS
     *
     * @param[in] nlocal number of orbitals
     * @param[in] Stmp overlap matrix
     * @param[in] Htmp H(t+dt/2) or H(t+dt)
     * @param[in] H_laststep H(t)
     * @param[in] print_matirx print internal matrix or not
     * @param[out] U_operator operator of propagator
     */
    void compute_propagator_etrs(const int nlocal,
                                 const std::complex<double>* Stmp,
                                 const std::complex<double>* Htmp,
                                 const std::complex<double>* H_laststep,
                                 std::complex<double>* U_operator,
                                 std::ofstream& ofs_running,
                                 const int print_matrix) const;
#endif // __MPI
};
} // namespace module_tddft

#endif
