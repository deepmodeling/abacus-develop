#include "edm_tddft.h"

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_lcao/module_rt/gather_mat.h"
#include "source_lcao/module_rt/propagator.h"

namespace module_dm
{

// Template function for EDM calculation supporting CPU and GPU
template <typename Device>
void edm_tddft_lapack(Parallel_Orbitals& pv,
                                 Setup_DM<std::complex<double>>& dmat,
                                 K_Vectors& kv,
                                 hamilt::Hamilt<std::complex<double>>* p_hamilt)
{
    ModuleBase::TITLE("elecstate", "edm_tddft_lapack");
    ModuleBase::timer::start("TD_Efficiency", "edm_tddft");

    // Dense gather/lapack/gemm operations below take the GLOBAL matrix
    // dimension; pv.nrow is only the per-process local row count.
    const int nlocal = pv.get_global_row_size();
    assert(nlocal >= 0);
    dmat.dm->edmk.resize(kv.get_nks());

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

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        p_hamilt->updateHk(ik);
        std::complex<double>* tmp_dmk_local = dmat.dm->get_dmk_ptr(ik);
        ModuleBase::ComplexMatrix& tmp_edmk = dmat.dm->edmk[ik];

#ifdef __MPI
        int myid = 0;
        const int root_proc = 0;
        int num_procs = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        // 1. Prepare Data Source Pointers (Host)
        // If np = 1, point directly to local data to avoid copy
        // If np > 1, gather data and point to the gathered buffer
        std::complex<double>* h_src = nullptr;
        std::complex<double>* s_src = nullptr;
        std::complex<double>* dmk_src = nullptr;

        // Global containers (Used only when num_procs > 1)
        module_rt::Matrix_g<std::complex<double>> h_mat_global, s_mat_global, dmk_global, edm_global;

        // Get Local Matrices
        hamilt::MatrixBlock<std::complex<double>> h_mat_local, s_mat_local;
        p_hamilt->matrix(h_mat_local, s_mat_local);

        if (num_procs == 1)
        {
            // Optimization: Direct access for single process
            h_src = h_mat_local.p;
            s_src = s_mat_local.p;
            dmk_src = tmp_dmk_local;
        }
        else
        {
            // Standard Gather Logic for multi-process
            module_rt::gatherMatrix(myid, root_proc, h_mat_local, h_mat_global);
            module_rt::gatherMatrix(myid, root_proc, s_mat_local, s_mat_global);

            hamilt::MatrixBlock<std::complex<double>> dmk_local_block;
            dmk_local_block.p = tmp_dmk_local;
            dmk_local_block.desc = pv.desc;
            module_rt::gatherMatrix(myid, root_proc, dmk_local_block, dmk_global);

            if (myid == root_proc)
            {
                h_src = h_mat_global.p.get();
                s_src = s_mat_global.p.get();
                dmk_src = dmk_global.p.get();
            }
        }

        // 2. GPU Calculation (on Rank 0)
        if (myid == root_proc)
        {
            ct::Tensor H_dev, S_dev, DMK_dev, ipiv_dev;

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
            ct::Tensor Sinv_dev = module_rt::create_identity_matrix<std::complex<double>>(nlocal, ct_device_type);

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

            // 3. Retrieve Result (D2H)
            std::complex<double>* edm_dest = nullptr;

            if (num_procs == 1)
            {
                // Directly copy to target local matrix
                tmp_edmk.create(pv.ncol, pv.nrow);
                edm_dest = tmp_edmk.c;
            }
            else
            {
                // Wait to set up edm_dest after allocating global buffer
                if (myid == root_proc && edm_global.p == nullptr)
                {
                    edm_global.p.reset(new std::complex<double>[nlocal * nlocal]);
                }
                edm_dest = edm_global.p.get();
            }

            if (num_procs == 1 || myid == root_proc)
            {
                syncmem_complex_d2h_op()(edm_dest, tmp4_dev.template data<std::complex<double>>(), nlocal * nlocal);
            }
        }

        // 4. Distribute (Only needed if num_procs > 1)
        if (num_procs > 1)
        {
            if (edm_global.p == nullptr)
            {
                edm_global.p.reset(new std::complex<double>[nlocal * nlocal]);
            }

            edm_global.row = nlocal;
            edm_global.col = nlocal;
            edm_global.desc.reset(new int[9]{1, pv.desc[1], nlocal, nlocal, nlocal, nlocal, 0, 0, nlocal});

            tmp_edmk.create(pv.ncol, pv.nrow);
            hamilt::MatrixBlock<std::complex<double>> edm_local_block;
            edm_local_block.p = tmp_edmk.c;
            edm_local_block.desc = pv.desc;
            module_rt::distributeMatrix(edm_local_block, edm_global);
        }
#else
        ModuleBase::WARNING_QUIT("elecstate::edm_tddft_lapack", "MPI is required for this function!");
#endif // __MPI
    } // end ik

#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Destroy cuBLAS & cuSOLVER handle
        ct::kernels::destroyGpuSolverHandle();
        ct::kernels::destroyGpuBlasHandle();
    }
#endif // __CUDA

    ModuleBase::timer::end("TD_Efficiency", "edm_tddft");
    return;
} // edm_tddft_lapack

// Explicit instantiation of template functions
template void edm_tddft_lapack<base_device::DEVICE_CPU>(Parallel_Orbitals& pv,
                                                                   Setup_DM<std::complex<double>>& dmat,
                                                                   K_Vectors& kv,
                                                                   hamilt::Hamilt<std::complex<double>>* p_hamilt);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void edm_tddft_lapack<base_device::DEVICE_GPU>(Parallel_Orbitals& pv,
                                                                   Setup_DM<std::complex<double>>& dmat,
                                                                   K_Vectors& kv,
                                                                   hamilt::Hamilt<std::complex<double>>* p_hamilt);
#endif // __CUDA

} // namespace module_dm
