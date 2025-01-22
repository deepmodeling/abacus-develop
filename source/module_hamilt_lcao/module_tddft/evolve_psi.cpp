#include "evolve_psi.h"

#include "band_energy.h"
#include "middle_hamilt.h"
#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"   // cuBLAS handle
#include "module_base/module_container/ATen/kernels/lapack.h" // cuSOLVER handle
#include "module_base/scalapack_connector.h"
#include "module_esolver/esolver_ks_lcao_tddft.h" // use gatherMatrix
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "norm_psi.h"
#include "propagator.h"
#include "upsi.h"

#include <complex>

namespace module_tddft
{
void evolve_psi(const int nband,
                const int nlocal,
                const Parallel_Orbitals* pv,
                hamilt::Hamilt<std::complex<double>>* p_hamilt,
                std::complex<double>* psi_k,
                std::complex<double>* psi_k_laststep,
                std::complex<double>* H_laststep,
                std::complex<double>* S_laststep,
                double* ekb,
                int htype,
                int propagator,
                std::ofstream& ofs_running,
                const int print_matrix)
{
    ofs_running << " evolve_psi::start " << std::endl;

    ModuleBase::TITLE("Evolve_psi", "evolve_psi");
    time_t time_start = time(nullptr);
    ofs_running << " Start Time : " << ctime(&time_start);

#ifdef __MPI

    hamilt::MatrixBlock<std::complex<double>> h_mat, s_mat;
    p_hamilt->matrix(h_mat, s_mat);

    std::complex<double>* Stmp = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Stmp, pv->nloc);
    BlasConnector::copy(pv->nloc, s_mat.p, 1, Stmp, 1);

    std::complex<double>* Htmp = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Htmp, pv->nloc);
    BlasConnector::copy(pv->nloc, h_mat.p, 1, Htmp, 1);

    std::complex<double>* Hold = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Hold, pv->nloc);
    BlasConnector::copy(pv->nloc, h_mat.p, 1, Hold, 1);

    std::complex<double>* U_operator = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(U_operator, pv->nloc);

    // (1)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute H(t+dt/2)
    /// @input H_laststep, Htmp, print_matrix
    /// @output Htmp
    if (htype == 1 && propagator != 2)
    {
        half_Hmatrix(pv, nband, nlocal, Htmp, Stmp, H_laststep, S_laststep, ofs_running, print_matrix);
    }

    // (2)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute U_operator
    /// @input Stmp, Htmp, print_matrix
    /// @output U_operator
    Propagator prop(propagator, pv, PARAM.mdp.md_dt);
    prop.compute_propagator(nlocal, Stmp, Htmp, H_laststep, U_operator, ofs_running, print_matrix);

    // (3)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief apply U_operator to the wave function of the previous step for new wave function
    /// @input U_operator, psi_k_laststep, print_matrix
    /// @output psi_k
    upsi(pv, nband, nlocal, U_operator, psi_k_laststep, psi_k, ofs_running, print_matrix);

    // (4)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief normalize psi_k
    /// @input Stmp, psi_not_norm, psi_k, print_matrix
    /// @output psi_k
    norm_psi(pv, nband, nlocal, Stmp, psi_k, ofs_running, print_matrix);

    // (5)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute ekb
    /// @input Htmp, psi_k
    /// @output ekb
    compute_ekb(pv, nband, nlocal, Hold, psi_k, ekb, ofs_running);

    delete[] Stmp;
    delete[] Htmp;
    delete[] Hold;
    delete[] U_operator;

#endif

    time_t time_end = time(nullptr);
    ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);

    ofs_running << " evolve_psi::end " << std::endl;

    return;
}

template <typename Device>
void evolve_psi_tensor(const int nband,
                       const int nlocal,
                       const Parallel_Orbitals* pv,
                       hamilt::Hamilt<std::complex<double>>* p_hamilt,
                       ct::Tensor& psi_k,
                       ct::Tensor& psi_k_laststep,
                       ct::Tensor& H_laststep,
                       ct::Tensor& S_laststep,
                       ct::Tensor& ekb,
                       int htype,
                       int propagator,
                       std::ofstream& ofs_running,
                       const int print_matrix,
                       const bool use_lapack)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    // Memory operations
    using syncmem_complex_h2d_op
        = base_device::memory::synchronize_memory_op<std::complex<double>, Device, base_device::DEVICE_CPU>;

#if ((defined __CUDA) /* || (defined __ROCM) */)
    // Initialize cuBLAS & cuSOLVER handle
    ct::kernels::createGpuSolverHandle();
    ct::kernels::createGpuBlasHandle();
#endif // __CUDA

    ofs_running << " evolve_psi_tensor::start " << std::endl;

    ModuleBase::TITLE("Evolve_psi", "evolve_psi");
    time_t time_start = time(nullptr);
    ofs_running << " Start Time : " << ctime(&time_start);

#ifdef __MPI

    hamilt::MatrixBlock<std::complex<double>> h_mat, s_mat;
    p_hamilt->matrix(h_mat, s_mat);

    // Create Tensor objects for temporary data and sync from host to device
    const int len_HS = use_lapack ? nlocal * nlocal : pv->nloc;
    ct::Tensor Stmp(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));
    ct::Tensor Htmp(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));
    ct::Tensor Hold(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));

    if (use_lapack)
    {
        // Need to gather H and S matrix to root process here
        int myid = 0;
        int num_procs = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        ModuleESolver::Matrix_g<std::complex<double>> h_mat_g, s_mat_g; // Global matrix structure

        // Collect H matrix
        ModuleESolver::gatherMatrix(myid, 0, h_mat, h_mat_g);
        syncmem_complex_h2d_op()(Htmp.data<std::complex<double>>(), h_mat_g.p.get(), len_HS);
        syncmem_complex_h2d_op()(Hold.data<std::complex<double>>(), h_mat_g.p.get(), len_HS);

        // Collect S matrix
        ModuleESolver::gatherMatrix(myid, 0, s_mat, s_mat_g);
        syncmem_complex_h2d_op()(Stmp.data<std::complex<double>>(), s_mat_g.p.get(), len_HS);
    }
    else
    {
        // Original code
        syncmem_complex_h2d_op()(Stmp.data<std::complex<double>>(), s_mat.p, len_HS);
        syncmem_complex_h2d_op()(Htmp.data<std::complex<double>>(), h_mat.p, len_HS);
        syncmem_complex_h2d_op()(Hold.data<std::complex<double>>(), h_mat.p, len_HS);
    }

    ct::Tensor U_operator(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({len_HS}));
    U_operator.zero();

    int myid = 0;
    int root_proc = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // (1)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute H(t+dt/2)
    /// @input H_laststep, Htmp, print_matrix
    /// @output Htmp
    if (htype == 1 && propagator != 2)
    {
        if (!use_lapack)
        {
            half_Hmatrix_tensor(pv, nband, nlocal, Htmp, Stmp, H_laststep, S_laststep, ofs_running, print_matrix);
        }
        else
        {
            if (myid == root_proc)
            {
                half_Hmatrix_tensor_lapack<Device>(pv,
                                                   nband,
                                                   nlocal,
                                                   Htmp,
                                                   Stmp,
                                                   H_laststep,
                                                   S_laststep,
                                                   ofs_running,
                                                   print_matrix);
            }
        }
    }

    // (2)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute U_operator
    /// @input Stmp, Htmp, print_matrix
    /// @output U_operator
    Propagator prop(propagator, pv, PARAM.mdp.md_dt);
    prop.compute_propagator_tensor<Device>(nlocal,
                                           Stmp,
                                           Htmp,
                                           H_laststep,
                                           U_operator,
                                           ofs_running,
                                           print_matrix,
                                           use_lapack);

    // (3)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief apply U_operator to the wave function of the previous step for new wave function
    /// @input U_operator, psi_k_laststep, print_matrix
    /// @output psi_k
    if (!use_lapack)
    {
        upsi_tensor(pv, nband, nlocal, U_operator, psi_k_laststep, psi_k, ofs_running, print_matrix);
    }
    else
    {
        if (myid == root_proc)
        {
            upsi_tensor_lapack<Device>(pv, nband, nlocal, U_operator, psi_k_laststep, psi_k, ofs_running, print_matrix);
        }
    }

    // (4)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief normalize psi_k
    /// @input Stmp, psi_not_norm, psi_k, print_matrix
    /// @output psi_k
    if (!use_lapack)
    {
        norm_psi_tensor(pv, nband, nlocal, Stmp, psi_k, ofs_running, print_matrix);
    }
    else
    {
        if (myid == root_proc)
        {
            norm_psi_tensor_lapack<Device>(pv, nband, nlocal, Stmp, psi_k, ofs_running, print_matrix);
        }
    }

    // (5)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute ekb
    /// @input Htmp, psi_k
    /// @output ekb
    if (!use_lapack)
    {
        compute_ekb_tensor(pv, nband, nlocal, Hold, psi_k, ekb, ofs_running);
    }
    else
    {
        if (myid == root_proc)
        {
            compute_ekb_tensor_lapack<Device>(pv, nband, nlocal, Hold, psi_k, ekb, ofs_running);
        }
    }

#endif // __MPI

    time_t time_end = time(nullptr);
    ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);

    ofs_running << " evolve_psi_tensor::end " << std::endl;

#if ((defined __CUDA) /* || (defined __ROCM) */)
    // Destroy cuBLAS & cuSOLVER handle
    ct::kernels::destroyGpuSolverHandle();
    ct::kernels::destroyGpuBlasHandle();
#endif // __CUDA

    return;
}

// Explicit instantiation of template functions
template void evolve_psi_tensor<base_device::DEVICE_CPU>(const int nband,
                                                         const int nlocal,
                                                         const Parallel_Orbitals* pv,
                                                         hamilt::Hamilt<std::complex<double>>* p_hamilt,
                                                         ct::Tensor& psi_k,
                                                         ct::Tensor& psi_k_laststep,
                                                         ct::Tensor& H_laststep,
                                                         ct::Tensor& S_laststep,
                                                         ct::Tensor& ekb,
                                                         int htype,
                                                         int propagator,
                                                         std::ofstream& ofs_running,
                                                         const int print_matrix,
                                                         const bool use_lapack);

#if ((defined __CUDA) /* || (defined __ROCM) */)
template void evolve_psi_tensor<base_device::DEVICE_GPU>(const int nband,
                                                         const int nlocal,
                                                         const Parallel_Orbitals* pv,
                                                         hamilt::Hamilt<std::complex<double>>* p_hamilt,
                                                         ct::Tensor& psi_k,
                                                         ct::Tensor& psi_k_laststep,
                                                         ct::Tensor& H_laststep,
                                                         ct::Tensor& S_laststep,
                                                         ct::Tensor& ekb,
                                                         int htype,
                                                         int propagator,
                                                         std::ofstream& ofs_running,
                                                         const int print_matrix,
                                                         const bool use_lapack);
#endif // __CUDA

} // namespace module_tddft