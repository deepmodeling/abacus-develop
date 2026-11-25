#include "edm_tddft.h"

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

namespace module_dm
{

// Template function for EDM calculation supporting CPU and GPU
template <typename Device>
void edm_tddft_lapack(const Parallel_Blkcyclic& layout,
                      const ModuleBase::MatrixBlock<std::complex<double>>& h_mat,
                      const ModuleBase::MatrixBlock<std::complex<double>>& s_mat,
                      std::complex<double>* dmk,
                      ModuleBase::ComplexMatrix& edmk)
{
    ModuleBase::TITLE("elecstate", "edm_tddft_lapack");
    ModuleBase::timer::start("TD_Efficiency", "edm_tddft_lapack");

    // Dense gather/lapack/gemm operations below take the GLOBAL matrix
    // dimension rather than the per-process local row count.
    const int nlocal = layout.nrglb();
    assert(nlocal >= 0);

    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Memory operations
    using syncmem_complex_h2d_op
        = base_device::memory::synchronize_memory_op<std::complex<double>, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op
        = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, Device>;

#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Initialize cuBLAS & cuSOLVER handle
        ct::kernels::createGpuSolverHandle();
        ct::kernels::createGpuBlasHandle();
    }
#endif // __CUDA

#ifdef __MPI
    const int nprocs = layout.nrproc() * layout.ncproc();
    const bool is_root = layout.blacs_in_this_processor(0, 0);
    const int matrix_size = nlocal * nlocal;
    const int context = layout.ictxt();
    // A single block on BLACS (0, 0), matching the original gather/scatter layout.
    int dense_desc[9] = {1, context, nlocal, nlocal, nlocal, nlocal, 0, 0, nlocal};
    std::array<int, 9> local_desc;
    std::copy_n(layout.desc(), 9, local_desc.begin());

    const std::complex<double>* h_src = h_mat.p;
    const std::complex<double>* s_src = s_mat.p;
    const std::complex<double>* dmk_src = dmk;
    std::vector<std::complex<double>> h_global;
    std::vector<std::complex<double>> s_global;
    std::vector<std::complex<double>> dmk_global;
    std::vector<std::complex<double>> edm_global;

    if (nprocs > 1)
    {
        h_global.resize(matrix_size);
        s_global.resize(matrix_size);
        dmk_global.resize(matrix_size);
        edm_global.resize(matrix_size);
        std::array<int, 9> h_desc;
        std::array<int, 9> s_desc;
        std::copy_n(h_mat.desc, 9, h_desc.begin());
        std::copy_n(s_mat.desc, 9, s_desc.begin());
        Cpxgemr2d(nlocal, nlocal, h_mat.p, 1, 1, h_desc.data(), h_global.data(), 1, 1, dense_desc, context);
        Cpxgemr2d(nlocal, nlocal, s_mat.p, 1, 1, s_desc.data(), s_global.data(), 1, 1, dense_desc, context);
        Cpxgemr2d(nlocal, nlocal, dmk, 1, 1, local_desc.data(), dmk_global.data(), 1, 1, dense_desc, context);
        h_src = h_global.data();
        s_src = s_global.data();
        dmk_src = dmk_global.data();
    }

    edmk.create(layout.ncloc(), layout.nrloc());
    // Compute on the BLACS root, on the selected device.
    if (is_root)
    {
        ct::Tensor H_dev;
        ct::Tensor S_dev;
        ct::Tensor DMK_dev;
        ct::Tensor ipiv_dev;

        // Allocate and Copy (H2D)
        H_dev = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        syncmem_complex_h2d_op()(H_dev.template data<std::complex<double>>(), h_src, nlocal * nlocal);

        S_dev = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        syncmem_complex_h2d_op()(S_dev.template data<std::complex<double>>(), s_src, nlocal * nlocal);

        DMK_dev = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        syncmem_complex_h2d_op()(DMK_dev.template data<std::complex<double>>(), dmk_src, nlocal * nlocal);

        ipiv_dev = ct::Tensor(ct::DataType::DT_INT, ct_device_type, ct::TensorShape({nlocal}));
        ipiv_dev.zero();

        // --- Calculate S^-1 using getrf + getrs ---
        // 1. LU decomposition S = P * L * U
        ct::kernels::lapack_getrf<std::complex<double>, ct_Device>()(nlocal,
                                                                     nlocal,
                                                                     S_dev.template data<std::complex<double>>(),
                                                                     nlocal,
                                                                     ipiv_dev.template data<int>());

        // 2. Solve S * Sinv = I
        ct::Tensor Sinv_dev(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        Sinv_dev.zero();
        const std::complex<double> diagonal_value = {1.0, 0.0};
        for (int i = 0; i < nlocal; ++i)
        {
            std::complex<double>* diagonal = Sinv_dev.template data<std::complex<double>>() + i * nlocal + i;
            ct::kernels::set_memory<std::complex<double>, ct_Device>()(diagonal, diagonal_value, 1);
        }

        ct::kernels::lapack_getrs<std::complex<double>, ct_Device>()('N',
                                                                     nlocal,
                                                                     nlocal,
                                                                     S_dev.template data<std::complex<double>>(),
                                                                     nlocal,
                                                                     ipiv_dev.template data<int>(),
                                                                     Sinv_dev.template data<std::complex<double>>(),
                                                                     nlocal);

        // --- EDM Calculation ---
        std::complex<double> one = {1.0, 0.0};
        std::complex<double> zero = {0.0, 0.0};

        // tmp1 = H * Sinv
        ct::Tensor tmp1_dev(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                                  'N',
                                                                  nlocal,
                                                                  nlocal,
                                                                  nlocal,
                                                                  &one,
                                                                  H_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  Sinv_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  &zero,
                                                                  tmp1_dev.template data<std::complex<double>>(),
                                                                  nlocal);

        // tmp2 = tmp1^T * DMK
        ct::Tensor tmp2_dev(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('T',
                                                                  'N',
                                                                  nlocal,
                                                                  nlocal,
                                                                  nlocal,
                                                                  &one,
                                                                  tmp1_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  DMK_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  &zero,
                                                                  tmp2_dev.template data<std::complex<double>>(),
                                                                  nlocal);

        // tmp3 = Sinv * H
        ct::Tensor tmp3_dev(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                                  'N',
                                                                  nlocal,
                                                                  nlocal,
                                                                  nlocal,
                                                                  &one,
                                                                  Sinv_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  H_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  &zero,
                                                                  tmp3_dev.template data<std::complex<double>>(),
                                                                  nlocal);

        // tmp4 = DMK * tmp3^T
        ct::Tensor tmp4_dev(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal, nlocal}));
        ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                                  'T',
                                                                  nlocal,
                                                                  nlocal,
                                                                  nlocal,
                                                                  &one,
                                                                  DMK_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  tmp3_dev.template data<std::complex<double>>(),
                                                                  nlocal,
                                                                  &zero,
                                                                  tmp4_dev.template data<std::complex<double>>(),
                                                                  nlocal);

        // tmp4 = tmp2 + tmp4
        ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                                  &one,
                                                                  tmp2_dev.template data<std::complex<double>>(),
                                                                  1,
                                                                  tmp4_dev.template data<std::complex<double>>(),
                                                                  1);

        // tmp4 = 0.5 * tmp4
        std::complex<double> half = {0.5, 0.0};
        ct::kernels::blas_scal<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                                  &half,
                                                                  tmp4_dev.template data<std::complex<double>>(),
                                                                  1);

        std::complex<double>* edm_dest = edmk.c;
        if (nprocs > 1)
        {
            edm_dest = edm_global.data();
        }
        syncmem_complex_d2h_op()(edm_dest, tmp4_dev.template data<std::complex<double>>(), matrix_size);
    }

    if (nprocs > 1)
    {
        Cpxgemr2d(nlocal, nlocal, edm_global.data(), 1, 1, dense_desc, edmk.c, 1, 1, local_desc.data(), context);
    }
#else
    ModuleBase::WARNING_QUIT("elecstate::edm_tddft_lapack", "MPI is required for this function!");
#endif // __MPI

#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Destroy cuBLAS & cuSOLVER handle
        ct::kernels::destroyGpuSolverHandle();
        ct::kernels::destroyGpuBlasHandle();
    }
#endif // __CUDA

    ModuleBase::timer::end("TD_Efficiency", "edm_tddft_lapack");
    return;
} // edm_tddft_lapack

// Explicit instantiation of template functions
template void edm_tddft_lapack<base_device::DEVICE_CPU>(const Parallel_Blkcyclic& layout,
    const ModuleBase::MatrixBlock<std::complex<double>>& h_mat,
    const ModuleBase::MatrixBlock<std::complex<double>>& s_mat,
    std::complex<double>* dmk,
    ModuleBase::ComplexMatrix& edmk);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void edm_tddft_lapack<base_device::DEVICE_GPU>(const Parallel_Blkcyclic& layout,
    const ModuleBase::MatrixBlock<std::complex<double>>& h_mat,
    const ModuleBase::MatrixBlock<std::complex<double>>& s_mat,
    std::complex<double>* dmk,
    ModuleBase::ComplexMatrix& edmk);
#endif // __CUDA

} // namespace module_dm
