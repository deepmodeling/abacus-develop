#include "cal_edm_tddft.h"

#include "source_base/module_container/ATen/core/tensor.h" // For ct::Tensor
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "source_base/module_device/memory_op.h"              // memory operations
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_io/module_parameter/parameter.h" // use PARAM.globalv

namespace elecstate
{
void print_local_matrix(std::ostream& os,
                        const std::complex<double>* matrix_data,
                        int local_rows, // pv.nrow
                        int local_cols, // pv.ncol
                        const std::string& matrix_name = "",
                        int rank = -1)
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
    ModuleBase::timer::tick("elecstate", "cal_edm_tddft");

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

    ModuleBase::timer::tick("elecstate", "cal_edm_tddft");
    return;
} // cal_edm_tddft

void cal_edm_tddft_tensor(Parallel_Orbitals& pv,
                          LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                          K_Vectors& kv,
                          hamilt::Hamilt<std::complex<double>>* p_hamilt)
{
    ModuleBase::timer::tick("elecstate", "cal_edm_tddft_tensor");

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

        // --- ScaLAPACK Inversion of S ---
        ct::Tensor ipiv_tensor(ct::DataType::DT_INT,
                               ct::DeviceType::CpuDevice,
                               ct::TensorShape({pv.nrow + pv.nb})); // Size for ScaLAPACK pivot array
        ipiv_tensor.zero();
        int* ipiv_ptr = ipiv_tensor.data<int>();

        int info = 0;
        const int one_int = 1;
        ScalapackConnector::getrf(nlocal, nlocal, Sinv_ptr, one_int, one_int, pv.desc, ipiv_ptr, &info);

        int lwork = -1;
        int liwork = -1;
        ct::Tensor work_query_tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({1}));
        ct::Tensor iwork_query_tensor(ct::DataType::DT_INT, ct::DeviceType::CpuDevice, ct::TensorShape({1}));

        ScalapackConnector::getri(nlocal,
                                  Sinv_ptr,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv_ptr,
                                  work_query_tensor.data<std::complex<double>>(),
                                  &lwork,
                                  iwork_query_tensor.data<int>(),
                                  &liwork,
                                  &info);

        // Resize work arrays based on query results
        lwork = work_query_tensor.data<std::complex<double>>()[0].real();
        work_query_tensor.resize(ct::TensorShape({lwork}));
        liwork = iwork_query_tensor.data<int>()[0];
        iwork_query_tensor.resize(ct::TensorShape({liwork}));

        ScalapackConnector::getri(nlocal,
                                  Sinv_ptr,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv_ptr,
                                  work_query_tensor.data<std::complex<double>>(),
                                  &lwork,
                                  iwork_query_tensor.data<int>(),
                                  &liwork,
                                  &info);

        // --- EDM Calculation using ScaLAPACK ---
        const char N_char = 'N';
        const char T_char = 'T';
        const std::complex<double> one_complex = {1.0, 0.0};
        const std::complex<double> zero_complex = {0.0, 0.0};
        const std::complex<double> half_complex = {0.5, 0.0};

        // tmp1 = Sinv * Htmp (result stored in tmp1)
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
                                 tmp1_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp2 = tmp1 * tmp_dmk (result stored in tmp2)
        ScalapackConnector::gemm(N_char,
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

        // tmp3 = Htmp * Sinv (result stored in tmp3)
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
                                 tmp3_ptr,
                                 one_int,
                                 one_int,
                                 pv.desc);

        // tmp4 = tmp_dmk * tmp3 (result stored in tmp4)
        ScalapackConnector::gemm(N_char,
                                 N_char,
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

        // tmp4 = 0.5 * tmp2 + 0.5 * tmp4 (final EDM contribution)
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
        // Serial version remains unchanged, using ModuleBase::ComplexMatrix directly
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
    ModuleBase::timer::tick("elecstate", "cal_edm_tddft_tensor");
    return;
} // cal_edm_tddft_tensor

} // namespace elecstate
