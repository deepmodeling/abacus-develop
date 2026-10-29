#ifdef __MPI
// #include <ATen/core/tensor_types.h>
#include "source_base/module_container/ATen/core/tensor_types.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "dm_diff.h"
namespace LR
{
    void CvX(
        const double* C,
        const BlockCyclicDistributionInformation& pc,
        const double* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CvX,
        const BlockCyclicDistributionInformation& pcx)
    {
        const int i1 = 1;
        const int ivirt = nocc + 1;
        char transa = 'N';
        char transb = 'N';
        const double alpha = 1.0;
        const double beta = 0;
        ScalapackConnector::gemm(transa,
            transb,
            naos,
            nocc,
            nvirt,
            alpha,
            C,
            i1,
            ivirt,
            pc.get_desc(),
            X,
            i1,
            i1,
            px.get_desc(),
            beta,
            CvX,
            i1,
            i1,
            pcx.get_desc());
    }
    void CvX(
        const std::complex<double>* C,
        const BlockCyclicDistributionInformation& pc,
        const std::complex<double>* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CvX,
        const BlockCyclicDistributionInformation& pcx)
    {
        const int i1 = 1;
        const int ivirt = nocc + 1;
        char transa = 'N';
        char transb = 'N';
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        ScalapackConnector::gemm(transa,
            transb,
            naos,
            nocc,
            nvirt,
            alpha,
            C,
            i1,
            ivirt,
            pc.get_desc(),
            X,
            i1,
            i1,
            px.get_desc(),
            beta,
            CvX,
            i1,
            i1,
            pcx.get_desc());
    }

    void CoXT(
        const double* C,
        const BlockCyclicDistributionInformation& pc,
        const double* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CoXT,
        const BlockCyclicDistributionInformation& pcxt)
    {
        const int i1 = 1;
        char transa = 'N';
        char transb = 'T';
        const double alpha = 1.0;
        const double beta = 0;
        ScalapackConnector::gemm(transa,
            transb,
            naos,
            nvirt,
            nocc,
            alpha,
            C,
            i1,
            i1,
            pc.get_desc(),
            X,
            i1,
            i1,
            px.get_desc(),
            beta,
            CoXT,
            i1,
            i1,
            pcxt.get_desc());
    }
    void CoXT(
        const std::complex<double>* C,
        const BlockCyclicDistributionInformation& pc,
        const std::complex<double>* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CoXT,
        const BlockCyclicDistributionInformation& pcxt)
    {
        const int i1 = 1;
        char transa = 'N';
        char transb = 'T';
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        ScalapackConnector::gemm(transa,
            transb,
            naos,
            nvirt,
            nocc,
            alpha,
            C,
            i1,
            i1,
            pc.get_desc(),
            X,
            i1,
            i1,
            px.get_desc(),
            beta,
            CoXT,
            i1,
            i1,
            pcxt.get_desc());
    }

    void AAT(
        const double* A,
        const BlockCyclicDistributionInformation& pa,
        const int& nrow,
        const int& ncol,
        double* C,
        const BlockCyclicDistributionInformation& pc,
        const bool add = false,
        const double& alpha = 1.0)
    {
        const int i1 = 1;
        char transa = 'N';
        char transb = 'T';
        const double beta = add ? 1.0 : 0.0;
        ScalapackConnector::gemm(transa,
            transb,
            nrow,
            nrow,
            ncol,
            alpha,
            A,
            i1,
            i1,
            pa.get_desc(),
            A,
            i1,
            i1,
            pa.get_desc(),
            beta,
            C,
            i1,
            i1,
            pc.get_desc());
    }
    void AAT(
        const std::complex<double>* A,
        const BlockCyclicDistributionInformation& pa,
        const int& nrow,
        const int& ncol,
        std::complex<double>* C,
        const BlockCyclicDistributionInformation& pc,
        const bool add = false,
        const std::complex<double>& alpha = std::complex<double>(1.0, 0.0))
    {
        const int i1 = 1;
        char transa = 'N';
        char transb = 'T';
        const std::complex<double> beta(add ? 1.0 : 0.0, 0.0);
        // take conjugate of A
        std::vector<std::complex<double>> A_conj(pa.get_local_size());
        for (int i = 0; i < pa.get_local_size(); ++i)
            A_conj[i] = std::conj(A[i]);
        ScalapackConnector::gemm(transa,
            transb,
            nrow,
            nrow,
            ncol,
            alpha,
            A_conj.data(),
            i1,
            i1,
            pa.get_desc(),
            A,
            i1,
            i1,
            pa.get_desc(),
            beta,
            C,
            i1,
            i1,
            pc.get_desc());
    }

    //output: col first, consistent with blas
    // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
    // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
    template<typename T>
    std::vector<ct::Tensor> cal_dm_diff_pblas(
        const T* const X_istate,
        const BlockCyclicDistributionInformation& px,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const BlockCyclicDistributionInformation& pc,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const BlockCyclicDistributionInformation& pmat,
        const bool renorm_k,
        const int nspin)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_diff_pblas");
        assert(px.comm() == pc.comm() && px.comm() == pmat.comm());
        assert(px.get_blacs_ctxt() == pc.get_blacs_ctxt() && px.get_blacs_ctxt() == pmat.get_blacs_ctxt());
        T empty_block(0);
        const T* const x_data = X_istate == nullptr ? &empty_block : X_istate;
        const int nks = c.get_nk();
        const int nk = nks / nspin;

        BlockCyclicDistributionInformation pcx;
        LR_Util::setup_2d_division(pcx, px.get_block_size(), naos, nocc, px.get_blacs_ctxt());
        ct::Tensor cvx(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { pcx.get_col_size(), pcx.get_row_size() });
        BlockCyclicDistributionInformation pcxt;
        LR_Util::setup_2d_division(pcxt, px.get_block_size(), naos, nvirt, px.get_blacs_ctxt());
        ct::Tensor coxt(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { pcxt.get_col_size(), pcxt.get_row_size() });
        std::vector<ct::Tensor> dm_diff(nks, ct::Tensor(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { pmat.get_col_size(), pmat.get_row_size() }));
        for (int iks = 0;iks < nks;++iks)
        {
            c.fix_k(iks);
            const int start = iks * px.get_local_size();
            // 1. C_virt * X
            CvX(c.get_pointer(), pc, x_data + start, px, naos, nocc, nvirt, cvx.data<T>(), pcx);
            // 2. C_occ * X^T
            CoXT(c.get_pointer(), pc, x_data + start, px, naos, nocc, nvirt, coxt.data<T>(), pcxt);
            // print_colfirst(c.get_pointer(), "c_pblas", naos, nocc + nvirt);
            // print_colfirst(X_istate.get_pointer(), "X_pblas", nvirt, nocc);
            // print_colfirst(cvx.data<T>(), "cvx_pblas", naos, nocc);
            // print_colfirst(coxt.data<T>(), "coxt_pblas", naos, nvirt);
            // 3. cvx*cvx^T - coxt*coxt^T
            AAT(cvx.data<T>(), pcx, naos, nocc, dm_diff[iks].data<T>(), pmat, false, renorm_k ? (T)(1.0 / (double)nk) : (T)1.0);
            // print_colfirst(dm_diff[iks].data<T>(), "dm_diff_1_pblas", naos, naos);
            AAT(coxt.data<T>(), pcxt, naos, nvirt, dm_diff[iks].data<T>(), pmat, true, renorm_k ? (T)(-1.0 / (double)nk) : (T)(-1.0));
            // print_colfirst(dm_diff[iks].data<T>(), "dm_diff_2_pblas", naos, naos);
        }
        return dm_diff;
    }
template std::vector<container::Tensor> cal_dm_diff_pblas<double>(const double* const, const BlockCyclicDistributionInformation&, const psi::Psi<double, base_device::DEVICE_CPU>&, const BlockCyclicDistributionInformation&, const int&, const int&, const int&, const BlockCyclicDistributionInformation&, const bool, const int);
template std::vector<container::Tensor> cal_dm_diff_pblas<std::complex<double>>(const std::complex<double>* const, const BlockCyclicDistributionInformation&, const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>&, const BlockCyclicDistributionInformation&, const int&, const int&, const int&, const BlockCyclicDistributionInformation&, const bool, const int);

}
#endif
