#include "evolve_psi.h"

#include "band_energy.h"
#include "gather_mat.h"
#include "middle_hamilt.h"
#include "norm_psi.h"
#include "propagator.h"
#include "solve_propagation.h"
#include "source_base/global_function.h"
#include "source_base/module_container/ATen/kernels/blas.h"   // cuBLAS handle
#include "source_base/module_container/ATen/kernels/lapack.h" // cuSOLVER handle
#include "source_base/module_external/blas_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_io/module_parameter/parameter.h"
#include "upsi.h"

#include "source_base/parallel_blkcyc.h"

#include <complex>
#include <vector>

namespace module_rt
{
void evolve_psi(const int nband,
                const int nlocal,
                const Parallel_Blkcyclic& matrix_layout,
                const Parallel_Blkcyclic& wfc_layout,
                const Parallel_Blkcyclic& mo_repr_layout,
                hamilt::Hamilt<std::complex<double>>* p_hamilt,
                std::complex<double>* psi_k,
                std::complex<double>* psi_k_laststep,
                std::complex<double>* H_laststep,
                std::complex<double>* S_laststep,
                std::complex<double>* P_k,
                const bool use_td_moving_gauge,
                double* ekb,
                int propagator,
                std::ofstream& ofs_running,
                const int print_matrix,
                const EijOutputInfo& output_info)
{
    ModuleBase::TITLE("module_rt", "evolve_psi");
    time_t time_start = time(nullptr);

#ifdef __MPI

    hamilt::MatrixBlock<std::complex<double>> h_mat;
    hamilt::MatrixBlock<std::complex<double>> s_mat;
    p_hamilt->matrix(h_mat, s_mat);

    std::complex<double>* Stmp = new std::complex<double>[matrix_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(Stmp, matrix_layout.size());
    BlasConnector::copy(matrix_layout.size(), s_mat.p, 1, Stmp, 1);

    std::complex<double>* Htmp = new std::complex<double>[matrix_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(Htmp, matrix_layout.size());
    BlasConnector::copy(matrix_layout.size(), h_mat.p, 1, Htmp, 1);

    std::complex<double>* Hold = new std::complex<double>[matrix_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(Hold, matrix_layout.size());
    BlasConnector::copy(matrix_layout.size(), h_mat.p, 1, Hold, 1);

    std::vector<std::complex<double>> U_operator(matrix_layout.size());
    Matrix2dView<const std::complex<double>> overlap_view(
        Stmp,
        matrix_layout.nrloc(),
        matrix_layout.ncloc(),
        matrix_layout.desc());
    Matrix2dView<const std::complex<double>> hamiltonian_view(
        Htmp,
        matrix_layout.nrloc(),
        matrix_layout.ncloc(),
        matrix_layout.desc());
    Matrix2dView<std::complex<double>> operator_view(
        U_operator.data(),
        matrix_layout.nrloc(),
        matrix_layout.ncloc(),
        matrix_layout.desc());

    // (1)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute H(t+dt/2)
    /// @input H_laststep, Htmp, print_matrix
    /// @output Htmp
    if (propagator != 2)
    {
        half_Hmatrix(matrix_layout, nband, nlocal, Htmp, Stmp, H_laststep, S_laststep, ofs_running, print_matrix);
    }

    // (2)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (propagator != 3)
    {
        /// @brief compute U_operator
        /// @input Stmp, Htmp, print_matrix
        /// @output U_operator
        Matrix2dView<const std::complex<double>> last_hamiltonian_view(
            H_laststep,
            matrix_layout.nrloc(),
            matrix_layout.ncloc(),
            matrix_layout.desc());
        const Matrix2dView<const std::complex<double>>* last_hamiltonian
            = propagator == 2 ? &last_hamiltonian_view : nullptr;
        Propagator prop(propagator, PARAM.inp.td_dt);
        prop.compute_propagator(overlap_view,
                                hamiltonian_view,
                                last_hamiltonian,
                                operator_view,
                                ofs_running,
                                print_matrix);
    }

    // (3)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (propagator != 3)
    {
        /// @brief apply U_operator to the wave function of the previous step for new wave function
        /// @input U_operator, psi_k_laststep, print_matrix
        /// @output psi_k
        upsi(matrix_layout,
             wfc_layout,
             nband,
             nlocal,
             U_operator.data(),
             psi_k_laststep,
             psi_k,
             ofs_running,
             print_matrix);
    }
    else
    {
        /// @brief solve the propagation equation
        /// @input Stmp, Htmp, psi_k_laststep
        /// @output psi_k
        Matrix2dView<const std::complex<double>> psi_last_view(
            psi_k_laststep,
            wfc_layout.nrloc(),
            wfc_layout.ncloc(),
            wfc_layout.desc());
        Matrix2dView<std::complex<double>> psi_view(
            psi_k, wfc_layout.nrloc(), wfc_layout.ncloc(), wfc_layout.desc());
        if (use_td_moving_gauge)
        {
            Matrix2dView<const std::complex<double>> gauge_view(
                P_k,
                matrix_layout.nrloc(),
                matrix_layout.ncloc(),
                matrix_layout.desc());
            solve_propagation(PARAM.inp.td_dt,
                              overlap_view,
                              hamiltonian_view,
                              gauge_view,
                              psi_last_view,
                              psi_view);
        }
        else
        {
            solve_propagation(PARAM.inp.td_dt,
                              overlap_view,
                              hamiltonian_view,
                              psi_last_view,
                              psi_view);
        }
    }

    // (4)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief normalize psi_k
    /// @input S(t+dt), psi_not_norm, psi_k, print_matrix
    /// @output psi_k
    norm_psi(matrix_layout,
             wfc_layout,
             mo_repr_layout,
             nband,
             nlocal,
             s_mat.p,
             psi_k,
             ofs_running,
             print_matrix);

    // (5)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute ekb
    /// @input Htmp, psi_k
    /// @output ekb
    compute_ekb(matrix_layout,
                wfc_layout,
                mo_repr_layout,
                nband,
                nlocal,
                Hold,
                psi_k,
                ekb,
                ofs_running,
                output_info);

    delete[] Stmp;
    delete[] Htmp;
    delete[] Hold;

#endif // __MPI

    time_t time_end = time(nullptr);
    ModuleBase::GlobalFunc::OUT_TIME("evolve_psi", time_start, time_end);

    return;
}

template <typename Device>
void evolve_psi_tensor(const int nband,
                       const int nlocal,
                       const Parallel_Blkcyclic& matrix_layout,
                       const Parallel_Blkcyclic& wfc_layout,
                       const Parallel_Blkcyclic& mo_repr_layout,
                       hamilt::Hamilt<std::complex<double>>* p_hamilt,
                       ct::Tensor& psi_k,
                       ct::Tensor& psi_k_laststep,
                       ct::Tensor& H_laststep,
                       ct::Tensor& S_laststep,
                       ct::Tensor& ekb,
                       int propagator,
                       std::ofstream& ofs_running,
                       const int print_matrix,
                       const bool use_lapack,
                       CublasMpResources& cublas_res,
                       const EijOutputInfo& output_info)
{
    ModuleBase::TITLE("module_rt", "evolve_psi_tensor");
    time_t time_start = time(nullptr);

    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    // Memory operations
    using syncmem_complex_h2d_op = base_device::memory::synchronize_memory_op<std::complex<double>, Device, base_device::DEVICE_CPU>;

#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Initialize cuBLAS & cuSOLVER handle
        ct::kernels::createGpuSolverHandle();
        ct::kernels::createGpuBlasHandle();
    }
#endif // __CUDA

#ifdef __MPI
    hamilt::MatrixBlock<std::complex<double>> h_mat, s_mat;
    p_hamilt->matrix(h_mat, s_mat);

    int myid = 0;
    int num_procs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    const int root_proc = 0;

    std::complex<double>* h_src = nullptr;
    std::complex<double>* s_src = nullptr;

    module_rt::Matrix_g<std::complex<double>> h_mat_g, s_mat_g;

    if (use_lapack)
    {
        if (num_procs == 1)
        {
            h_src = h_mat.p;
            s_src = s_mat.p;
        }
        else
        {
            module_rt::gatherMatrix(myid, 0, h_mat, h_mat_g);
            module_rt::gatherMatrix(myid, 0, s_mat, s_mat_g);
            if (myid == root_proc)
            {
                h_src = h_mat_g.p.get();
                s_src = s_mat_g.p.get();
            }
        }
    }
    else
    {
        h_src = h_mat.p;
        s_src = s_mat.p;
    }

    const int len_HS = use_lapack ? nlocal * nlocal : matrix_layout.size();

    ct::Tensor Stmp(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));

    if (s_src != nullptr)
    {
        if (!use_lapack || myid == root_proc)
        {
            ModuleBase::timer::start("TD_Efficiency", "host_device_comm");
            syncmem_complex_h2d_op()(Stmp.data<std::complex<double>>(), s_src, len_HS);
            ModuleBase::timer::end("TD_Efficiency", "host_device_comm");
        }
    }

    ct::Tensor Htmp(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));

    if (h_src != nullptr)
    {
        if (!use_lapack || myid == root_proc)
        {
            ModuleBase::timer::start("TD_Efficiency", "host_device_comm");
            syncmem_complex_h2d_op()(Htmp.data<std::complex<double>>(), h_src, len_HS);
            ModuleBase::timer::end("TD_Efficiency", "host_device_comm");
        }
    }

    // (1) Compute H(t+dt/2)
    if (propagator != 2)
    {
        if (!use_lapack)
        {
            half_Hmatrix_tensor(matrix_layout, nband, nlocal, Htmp, Stmp, H_laststep, S_laststep, ofs_running, print_matrix, cublas_res);
        }
        else if (myid == root_proc)
        {
            half_Hmatrix_tensor_lapack<Device>(nband, nlocal, Htmp, Stmp, H_laststep, S_laststep, ofs_running, print_matrix);
        }
    }

    ct::Tensor U_operator(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));
    U_operator.zero();

    // (2) Compute U_operator
    Propagator prop(propagator, PARAM.inp.td_dt);
    prop.compute_propagator_tensor<Device>(matrix_layout,
                                           nlocal,
                                           Stmp,
                                           Htmp,
                                           H_laststep,
                                           U_operator,
                                           ofs_running,
                                           print_matrix,
                                           use_lapack,
                                           cublas_res);

    // (3) Apply U_operator (psi_k = U * psi_last)
    if (!use_lapack)
    {
        upsi_tensor(matrix_layout,
                    wfc_layout,
                    nband,
                    nlocal,
                    U_operator,
                    psi_k_laststep,
                    psi_k,
                    ofs_running,
                    print_matrix,
                    cublas_res);
    }
    else if (myid == root_proc)
    {
        upsi_tensor_lapack<Device>(nband, nlocal, U_operator, psi_k_laststep, psi_k, ofs_running, print_matrix);
    }

    // (4) Normalize in the endpoint metric, not the averaged propagation metric.
    if (s_src != nullptr && (!use_lapack || myid == root_proc))
    {
        syncmem_complex_h2d_op()(Stmp.data<std::complex<double>>(), s_src, len_HS);
    }
    if (!use_lapack)
    {
        norm_psi_tensor(matrix_layout,
                        wfc_layout,
                        mo_repr_layout,
                        nband,
                        nlocal,
                        Stmp,
                        psi_k,
                        ofs_running,
                        print_matrix,
                        cublas_res);
    }
    else if (myid == root_proc)
    {
        norm_psi_tensor_lapack<Device>(nband, nlocal, Stmp, psi_k, ofs_running, print_matrix);
    }

    // (5) Compute ekb
    ct::Tensor Hold(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));

    // Resync H matrix
    if (h_src != nullptr)
    {
        if (!use_lapack || myid == root_proc)
        {
            ModuleBase::timer::start("TD_Efficiency", "host_device_comm");
            syncmem_complex_h2d_op()(Hold.data<std::complex<double>>(), h_src, len_HS);
            ModuleBase::timer::end("TD_Efficiency", "host_device_comm");
        }
    }

    if (!use_lapack)
    {
        compute_ekb_tensor(matrix_layout,
                           wfc_layout,
                           mo_repr_layout,
                           nband,
                           nlocal,
                           Hold,
                           psi_k,
                           ekb,
                           ofs_running,
                           cublas_res);
    }
    else if (myid == root_proc)
    {
        compute_ekb_tensor_lapack<Device>(nband, nlocal, Hold, psi_k, ekb, ofs_running, output_info);
    }
#endif // __MPI

#if ((defined __CUDA) /* || (defined __ROCM) */)
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Destroy cuBLAS & cuSOLVER handle
        ct::kernels::destroyGpuSolverHandle();
        ct::kernels::destroyGpuBlasHandle();
    }
#endif // __CUDA

    time_t time_end = time(nullptr);
    ModuleBase::GlobalFunc::OUT_TIME("evolve_psi", time_start, time_end);

    return;
}

// Explicit instantiation of template functions
template void evolve_psi_tensor<base_device::DEVICE_CPU>(const int nband,
                                                         const int nlocal,
                                                         const Parallel_Blkcyclic& matrix_layout,
                                                         const Parallel_Blkcyclic& wfc_layout,
                                                         const Parallel_Blkcyclic& mo_repr_layout,
                                                         hamilt::Hamilt<std::complex<double>>* p_hamilt,
                                                         ct::Tensor& psi_k,
                                                         ct::Tensor& psi_k_laststep,
                                                         ct::Tensor& H_laststep,
                                                         ct::Tensor& S_laststep,
                                                         ct::Tensor& ekb,
                                                         int propagator,
                                                         std::ofstream& ofs_running,
                                                         const int print_matrix,
                                                         const bool use_lapack,
                                                         CublasMpResources& cublas_res,
                                                         const EijOutputInfo& output_info);

#if ((defined __CUDA) /* || (defined __ROCM) */)
template void evolve_psi_tensor<base_device::DEVICE_GPU>(const int nband,
                                                         const int nlocal,
                                                         const Parallel_Blkcyclic& matrix_layout,
                                                         const Parallel_Blkcyclic& wfc_layout,
                                                         const Parallel_Blkcyclic& mo_repr_layout,
                                                         hamilt::Hamilt<std::complex<double>>* p_hamilt,
                                                         ct::Tensor& psi_k,
                                                         ct::Tensor& psi_k_laststep,
                                                         ct::Tensor& H_laststep,
                                                         ct::Tensor& S_laststep,
                                                         ct::Tensor& ekb,
                                                         int propagator,
                                                         std::ofstream& ofs_running,
                                                         const int print_matrix,
                                                         const bool use_lapack,
                                                         CublasMpResources& cublas_res,
                                                         const EijOutputInfo& output_info);
#endif // __CUDA

} // namespace module_rt
