#include "cal_edm_tddft.h"

#include "source_base/module_container/ATen/core/tensor.h" // For ct::Tensor
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "source_base/module_device/memory_op.h"              // memory operations
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_io/module_parameter/parameter.h" // use PARAM.globalv
#include "source_lcao/module_rt/gather_mat.h"     // gatherMatrix and distributeMatrix
#include "source_lcao/module_rt/propagator.h"     // Include header for create_identity_matrix

namespace elecstate
{
void print_local_matrix(std::ostream& os,
                        const std::complex<double>* matrix_data,
                        int local_rows,
                        int local_cols,
                        const std::string& matrix_name,
                        int rank)
{
    if (!matrix_name.empty() || rank >= 0)
    {
        os << "=== ";
        if (!matrix_name.empty())
        {
            os << "Matrix: " << matrix_name;
            if (rank >= 0)
                os << " ";
        }
        if (rank >= 0)
        {
            os << "(Process: " << rank + 1 << ")";
        }
        os << " (Local dims: " << local_rows << " x " << local_cols << ") ===" << std::endl;
    }

    os << std::fixed << std::setprecision(10) << std::showpos;

    for (int i = 0; i < local_rows; ++i) // Iterate over rows (i)
    {
        for (int j = 0; j < local_cols; ++j) // Iterate over columns (j)
        {
            // For column-major storage, element (i, j) is at index i + j * LDA
            // where LDA (leading dimension) is typically the number of *rows* in the local block.
            int idx = i + j * local_rows;
            os << "(" << std::real(matrix_data[idx]) << "," << std::imag(matrix_data[idx]) << ") ";
        }
        os << std::endl; // New line after each row
    }
    os.unsetf(std::ios_base::fixed | std::ios_base::showpos);
    os << std::endl;
}

// use the original formula (Hamiltonian matrix) to calculate energy density matrix
void cal_edm_tddft(Parallel_Orbitals& pv,
                   LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                   K_Vectors& kv,
                   hamilt::Hamilt<std::complex<double>>* p_hamilt)
{
    ModuleBase::TITLE("elecstate", "cal_edm_tddft");
    ModuleBase::timer::tick("TD_Efficiency", "cal_edm_tddft");

    const int nlocal = PARAM.globalv.nlocal;
    assert(nlocal >= 0);

    dmat.dm->EDMK.resize(kv.get_nks());

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        p_hamilt->updateHk(ik);
        std::complex<double>* tmp_dmk = dmat.dm->get_DMK_pointer(ik);
        ModuleBase::ComplexMatrix& tmp_edmk = dmat.dm->EDMK[ik];

#ifdef __MPI
        const int nloc = pv.nloc;
        const int ncol = pv.ncol;
        const int nrow = pv.nrow;

        tmp_edmk.create(ncol, nrow);
        std::complex<double>* Htmp = new std::complex<double>[nloc];
        std::complex<double>* Sinv = new std::complex<double>[nloc];
        std::complex<double>* tmp1 = new std::complex<double>[nloc];
        std::complex<double>* tmp2 = new std::complex<double>[nloc];
        std::complex<double>* tmp3 = new std::complex<double>[nloc];
        std::complex<double>* tmp4 = new std::complex<double>[nloc];

        ModuleBase::GlobalFunc::ZEROS(Htmp, nloc);
        ModuleBase::GlobalFunc::ZEROS(Sinv, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp1, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp2, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp3, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp4, nloc);

        const int inc = 1;

        hamilt::MatrixBlock<std::complex<double>> h_mat;
        hamilt::MatrixBlock<std::complex<double>> s_mat;

        p_hamilt->matrix(h_mat, s_mat);
        BlasConnector::copy(nloc, h_mat.p, inc, Htmp, inc);
        BlasConnector::copy(nloc, s_mat.p, inc, Sinv, inc);

        vector<int> ipiv(nloc, 0);
        int info = 0;
        const int one_int = 1;

        ScalapackConnector::getrf(nlocal, nlocal, Sinv, one_int, one_int, pv.desc, ipiv.data(), &info);

        int lwork = -1;
        int liwork = -1;

        // if lwork == -1, then the size of work is (at least) of length 1.
        std::vector<std::complex<double>> work(1, 0);

        // if liwork = -1, then the size of iwork is (at least) of length 1.
        std::vector<int> iwork(1, 0);

        ScalapackConnector::getri(nlocal,
                                  Sinv,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv.data(),
                                  work.data(),
                                  &lwork,
                                  iwork.data(),
                                  &liwork,
                                  &info);

        lwork = work[0].real();
        work.resize(lwork, 0);
        liwork = iwork[0];
        iwork.resize(liwork, 0);

        ScalapackConnector::getri(nlocal,
                                  Sinv,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv.data(),
                                  work.data(),
                                  &lwork,
                                  iwork.data(),
                                  &liwork,
                                  &info);

        const char N_char = 'N';
        const char T_char = 'T';
        const std::complex<double> one_complex = {1.0, 0.0};
        const std::complex<double> zero_complex = {0.0, 0.0};
        const std::complex<double> half_complex = {0.5, 0.0};

        // tmp1 = Htmp * Sinv
        ScalapackConnector::gemm(N_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 Htmp,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 Sinv,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp1,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp2 = tmp1^T * tmp_dmk
        ScalapackConnector::gemm(T_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp1,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 tmp_dmk,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp2,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp3 = Sinv * Htmp
        ScalapackConnector::gemm(N_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 Sinv,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 Htmp,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp3,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp4 = tmp_dmk * tmp3^T
        ScalapackConnector::gemm(N_char,
                                 T_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp_dmk,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 tmp3,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp4,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp4 = 0.5 * (tmp2 + tmp4)
        ScalapackConnector::geadd(N_char,
                                  nlocal,
                                  nlocal,
                                  half_complex,
                                  tmp2,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  half_complex,
                                  tmp4,
                                  one_int,
                                  one_int,
                                  pv.desc);

        BlasConnector::copy(nloc, tmp4, inc, tmp_edmk.c, inc);

        delete[] Htmp;
        delete[] Sinv;
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] tmp4;
#else
        // for serial version
        tmp_edmk.create(pv.ncol, pv.nrow);
        ModuleBase::ComplexMatrix Sinv(nlocal, nlocal);
        ModuleBase::ComplexMatrix Htmp(nlocal, nlocal);

        hamilt::MatrixBlock<std::complex<double>> h_mat;
        hamilt::MatrixBlock<std::complex<double>> s_mat;

        p_hamilt->matrix(h_mat, s_mat);

        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                Htmp(i, j) = h_mat.p[i * nlocal + j];
                Sinv(i, j) = s_mat.p[i * nlocal + j];
            }
        }
        int INFO = 0;

        int lwork = 3 * nlocal - 1; // tmp
        std::complex<double>* work = new std::complex<double>[lwork];
        ModuleBase::GlobalFunc::ZEROS(work, lwork);

        int IPIV[nlocal];

        LapackConnector::zgetrf(nlocal, nlocal, Sinv, nlocal, IPIV, &INFO);
        LapackConnector::zgetri(nlocal, Sinv, nlocal, IPIV, work, lwork, &INFO);
        // I just use ModuleBase::ComplexMatrix temporarily, and will change it
        // to std::complex<double>*
        ModuleBase::ComplexMatrix tmp_dmk_base(nlocal, nlocal);
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                tmp_dmk_base(i, j) = tmp_dmk[i * nlocal + j];
            }
        }
        tmp_edmk = 0.5 * (Sinv * Htmp * tmp_dmk_base + tmp_dmk_base * Htmp * Sinv);
        delete[] work;
#endif
    } // end ik

    ModuleBase::timer::tick("TD_Efficiency", "cal_edm_tddft");
    return;
} // cal_edm_tddft

void cal_edm_tddft_tensor(Parallel_Orbitals& pv,
                          LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                          K_Vectors& kv,
                          hamilt::Hamilt<std::complex<double>>* p_hamilt)
{
    ModuleBase::TITLE("elecstate", "cal_edm_tddft_tensor");
    ModuleBase::timer::tick("TD_Efficiency", "cal_edm_tddft");

    const int nlocal = PARAM.globalv.nlocal;
    assert(nlocal >= 0);
    dmat.dm->EDMK.resize(kv.get_nks());

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        p_hamilt->updateHk(ik);
        std::complex<double>* tmp_dmk = dmat.dm->get_DMK_pointer(ik);
        ModuleBase::ComplexMatrix& tmp_edmk = dmat.dm->EDMK[ik];

#ifdef __MPI
        const int nloc = pv.nloc;
        const int ncol = pv.ncol;
        const int nrow = pv.nrow;

        // Initialize EDMK matrix
        tmp_edmk.create(ncol, nrow);

        // Allocate Tensor objects on CPU
        ct::Tensor Htmp_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nloc}));
        Htmp_tensor.zero();

        ct::Tensor Sinv_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nloc}));
        Sinv_tensor.zero();

        ct::Tensor tmp1_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nloc}));
        tmp1_tensor.zero();

        ct::Tensor tmp2_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nloc}));
        tmp2_tensor.zero();

        ct::Tensor tmp3_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nloc}));
        tmp3_tensor.zero();

        ct::Tensor tmp4_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nloc}));
        tmp4_tensor.zero();

        // Get raw pointers from tensors for ScaLAPACK calls
        std::complex<double>* Htmp_ptr = Htmp_tensor.data<std::complex<double>>();
        std::complex<double>* Sinv_ptr = Sinv_tensor.data<std::complex<double>>();
        std::complex<double>* tmp1_ptr = tmp1_tensor.data<std::complex<double>>();
        std::complex<double>* tmp2_ptr = tmp2_tensor.data<std::complex<double>>();
        std::complex<double>* tmp3_ptr = tmp3_tensor.data<std::complex<double>>();
        std::complex<double>* tmp4_ptr = tmp4_tensor.data<std::complex<double>>();

        const int inc = 1;
        hamilt::MatrixBlock<std::complex<double>> h_mat;
        hamilt::MatrixBlock<std::complex<double>> s_mat;
        p_hamilt->matrix(h_mat, s_mat);

        // Copy Hamiltonian and Overlap matrices into Tensor buffers using BlasConnector
        BlasConnector::copy(nloc, h_mat.p, inc, Htmp_ptr, inc);
        BlasConnector::copy(nloc, s_mat.p, inc, Sinv_ptr, inc);

        int myid = 0;
        const int root_proc = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        // --- ScaLAPACK Inversion of S ---
        ct::Tensor ipiv(ct::DataType::DT_INT,
                        ct::DeviceType::CpuDevice,
                        ct::TensorShape({pv.nrow + pv.nb})); // Size for ScaLAPACK pivot array
        ipiv.zero();
        int* ipiv_ptr = ipiv.data<int>();

        int info = 0;
        const int one_int = 1;
        ScalapackConnector::getrf(nlocal, nlocal, Sinv_ptr, one_int, one_int, pv.desc, ipiv_ptr, &info);

        int lwork = -1;
        int liwork = -1;
        ct::Tensor work_query(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({1}));
        ct::Tensor iwork_query(ct::DataType::DT_INT, ct::DeviceType::CpuDevice, ct::TensorShape({1}));

        ScalapackConnector::getri(nlocal,
                                  Sinv_ptr,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv_ptr,
                                  work_query.data<std::complex<double>>(),
                                  &lwork,
                                  iwork_query.data<int>(),
                                  &liwork,
                                  &info);

        // Resize work arrays based on query results
        lwork = work_query.data<std::complex<double>>()[0].real();
        work_query.resize(ct::TensorShape({lwork}));
        liwork = iwork_query.data<int>()[0];
        iwork_query.resize(ct::TensorShape({liwork}));

        ScalapackConnector::getri(nlocal,
                                  Sinv_ptr,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv_ptr,
                                  work_query.data<std::complex<double>>(),
                                  &lwork,
                                  iwork_query.data<int>(),
                                  &liwork,
                                  &info);

        // --- EDM Calculation using ScaLAPACK ---
        const char N_char = 'N';
        const char T_char = 'T';
        const std::complex<double> one_complex = {1.0, 0.0};
        const std::complex<double> zero_complex = {0.0, 0.0};
        const std::complex<double> half_complex = {0.5, 0.0};

        // tmp1 = Htmp * Sinv
        ScalapackConnector::gemm(N_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 Htmp_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 Sinv_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp1_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp2 = tmp1^T * tmp_dmk
        ScalapackConnector::gemm(T_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp1_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 tmp_dmk,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp2_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp3 = Sinv * Htmp
        ScalapackConnector::gemm(N_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 Sinv_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 Htmp_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp3_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp4 = tmp_dmk * tmp3^T
        ScalapackConnector::gemm(N_char,
                                 T_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp_dmk,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 tmp3_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp4_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp4 = 0.5 * (tmp2 + tmp4)
        ScalapackConnector::geadd(N_char,
                                  nlocal,
                                  nlocal,
                                  half_complex,
                                  tmp2_ptr,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  half_complex,
                                  tmp4_ptr,
                                  one_int,
                                  one_int,
                                  pv.desc);

        // Copy final result from Tensor buffer back to EDMK matrix
        BlasConnector::copy(nloc, tmp4_ptr, inc, tmp_edmk.c, inc);

#else
        ModuleBase::WARNING_QUIT("elecstate::cal_edm_tddft_tensor", "MPI is required for this function!");
#endif
    } // end ik
    ModuleBase::timer::tick("TD_Efficiency", "cal_edm_tddft");
    return;
} // cal_edm_tddft_tensor

// Template function for EDM calculation supporting CPU and GPU
template <typename Device>
void cal_edm_tddft_tensor_lapack(Parallel_Orbitals& pv,
                                 LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                                 K_Vectors& kv,
                                 hamilt::Hamilt<std::complex<double>>* p_hamilt)
{
    ModuleBase::TITLE("elecstate", "cal_edm_tddft_tensor_lapack");
    ModuleBase::timer::tick("TD_Efficiency", "cal_edm_tddft");

    const int nlocal = PARAM.globalv.nlocal;
    assert(nlocal >= 0);
    dmat.dm->EDMK.resize(kv.get_nks());

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
        std::complex<double>* tmp_dmk_local = dmat.dm->get_DMK_pointer(ik);
        ModuleBase::ComplexMatrix& tmp_edmk = dmat.dm->EDMK[ik];

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
            auto Sinv_dev = module_rt::create_identity_matrix<std::complex<double>>(nlocal, ct_device_type);

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
        ModuleBase::WARNING_QUIT("elecstate::cal_edm_tddft_tensor_lapack", "MPI is required for this function!");
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

    ModuleBase::timer::tick("TD_Efficiency", "cal_edm_tddft");
    return;
} // cal_edm_tddft_tensor_lapack

// Explicit instantiation of template functions
template void cal_edm_tddft_tensor_lapack<base_device::DEVICE_CPU>(Parallel_Orbitals& pv,
                                                                   LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                                                                   K_Vectors& kv,
                                                                   hamilt::Hamilt<std::complex<double>>* p_hamilt);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void cal_edm_tddft_tensor_lapack<base_device::DEVICE_GPU>(Parallel_Orbitals& pv,
                                                                   LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                                                                   K_Vectors& kv,
                                                                   hamilt::Hamilt<std::complex<double>>* p_hamilt);
#endif // __CUDA

} // namespace elecstate
