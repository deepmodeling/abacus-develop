#include "edm_tddft.h"

#include "source_base/module_container/ATen/core/tensor.h" // For ct::Tensor
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "source_base/module_device/memory_op.h"              // memory operations
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_lcao/module_rt/gather_mat.h"     // gatherMatrix and distributeMatrix
#include "source_lcao/module_rt/propagator.h"     // Include header for create_identity_matrix

namespace module_dm
{
// use the original formula (Hamiltonian matrix) to calculate energy density matrix
void edm_tddft(Parallel_Orbitals& pv,
                   Setup_DM<std::complex<double>>& dmat,
                   K_Vectors& kv,
                   hamilt::Hamilt<std::complex<double>>* p_hamilt)
{
    ModuleBase::TITLE("elecstate", "edm_tddft");
    ModuleBase::timer::start("TD_Efficiency", "edm_tddft");

    // Scalapack getrf/getri/gemm/geadd below take the GLOBAL matrix
    // dimension; pv.nrow is only the per-process local row count and
    // would truncate the operations to a submatrix of the full matrix.
    const int nlocal = pv.get_global_row_size();
    assert(nlocal >= 0);

    dmat.dm->edmk.resize(kv.get_nks());

    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        p_hamilt->updateHk(ik);
        std::complex<double>* tmp_dmk = dmat.dm->get_dmk_ptr(ik);
        ModuleBase::ComplexMatrix& tmp_edmk = dmat.dm->edmk[ik];

#ifdef __MPI
        const int nloc = pv.nloc;
        const int ncol = pv.ncol;
        const int nrow = pv.nrow;

        tmp_edmk.create(ncol, nrow);
        std::vector<std::complex<double>> Htmp_vec(nloc);
        std::vector<std::complex<double>> Sinv_vec(nloc);
        std::vector<std::complex<double>> tmp1_vec(nloc);
        std::vector<std::complex<double>> tmp2_vec(nloc);
        std::vector<std::complex<double>> tmp3_vec(nloc);
        std::vector<std::complex<double>> tmp4_vec(nloc);
        std::complex<double>* Htmp = Htmp_vec.data();
        std::complex<double>* Sinv = Sinv_vec.data();
        std::complex<double>* tmp1 = tmp1_vec.data();
        std::complex<double>* tmp2 = tmp2_vec.data();
        std::complex<double>* tmp3 = tmp3_vec.data();
        std::complex<double>* tmp4 = tmp4_vec.data();

        // The buffers are value-initialized to zero by the std::vector
        // constructors and fully overwritten below (Htmp/Sinv by copy,
        // tmp1~tmp4 by gemm with beta == 0), so no explicit zeroing is needed.
        const int inc = 1;

        hamilt::MatrixBlock<std::complex<double>> h_mat;
        hamilt::MatrixBlock<std::complex<double>> s_mat;

        p_hamilt->matrix(h_mat, s_mat);
        BlasConnector::copy(nloc, h_mat.p, inc, Htmp, inc);
        BlasConnector::copy(nloc, s_mat.p, inc, Sinv, inc);

        std::vector<int> ipiv(nloc, 0);
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
        std::vector<std::complex<double>> work_vec(lwork);
        std::complex<double>* work = work_vec.data();

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
#endif
    } // end ik

    ModuleBase::timer::end("TD_Efficiency", "edm_tddft");
    return;
} // edm_tddft

} // namespace module_dm
