#include "propagator.h"

#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_container/ATen/kernels/lapack.h"
#include "module_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "module_base/module_device/memory_op.h"              // memory operations
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
Propagator::~Propagator()
{
}
#ifdef __MPI
void Propagator::compute_propagator(const int nlocal,
                                    const std::complex<double>* Stmp,
                                    const std::complex<double>* Htmp,
                                    const std::complex<double>* H_laststep,
                                    std::complex<double>* U_operator,
                                    std::ofstream& ofs_running,
                                    const int print_matrix) const
{
    int tag;
    switch (ptype)
    {
    case 0:
        compute_propagator_cn2(nlocal, Stmp, Htmp, U_operator, ofs_running, print_matrix);
        break;

    case 1:
        tag = 1;
        compute_propagator_taylor(nlocal, Stmp, Htmp, U_operator, ofs_running, print_matrix, tag);
        break;

    case 2:
        compute_propagator_etrs(nlocal, Stmp, Htmp, H_laststep, U_operator, ofs_running, print_matrix);
        break;

    default:
        ModuleBase::WARNING_QUIT("Propagator::compute_propagator", "Method of RT-TDDFT propagator is wrong!");
        break;
    }
}

template <typename Device>
void Propagator::compute_propagator_tensor(const int nlocal,
                                           const ct::Tensor& Stmp,
                                           const ct::Tensor& Htmp,
                                           const ct::Tensor& H_laststep,
                                           ct::Tensor& U_operator,
                                           std::ofstream& ofs_running,
                                           const int print_matrix,
                                           const bool use_lapack) const
{
    int tag;
    switch (ptype)
    {
    case 0:
        if (!use_lapack)
        {
            compute_propagator_cn2_tensor(nlocal, Stmp, Htmp, U_operator, ofs_running, print_matrix);
        }
        else
        {
            int myid = 0;
            int root_proc = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &myid);
            if (myid == root_proc)
            {
                compute_propagator_cn2_tensor_lapack<Device>(nlocal, Stmp, Htmp, U_operator, ofs_running, print_matrix);
            }
        }
        break;

    default:
        ModuleBase::WARNING_QUIT("Propagator::compute_propagator_tensor",
                                 "The Tensor-based RT-TDDFT propagator currently supports Crank–Nicolson method only!");
        break;
    }
}

// Explicit instantiation of template functions
template void Propagator::compute_propagator_tensor<base_device::DEVICE_CPU>(const int nlocal,
                                                                             const ct::Tensor& Stmp,
                                                                             const ct::Tensor& Htmp,
                                                                             const ct::Tensor& H_laststep,
                                                                             ct::Tensor& U_operator,
                                                                             std::ofstream& ofs_running,
                                                                             const int print_matrix,
                                                                             const bool use_lapack) const;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void Propagator::compute_propagator_tensor<base_device::DEVICE_GPU>(const int nlocal,
                                                                             const ct::Tensor& Stmp,
                                                                             const ct::Tensor& Htmp,
                                                                             const ct::Tensor& H_laststep,
                                                                             ct::Tensor& U_operator,
                                                                             std::ofstream& ofs_running,
                                                                             const int print_matrix,
                                                                             const bool use_lapack) const;
#endif // __CUDA
#endif // __MPI
} // namespace module_tddft
