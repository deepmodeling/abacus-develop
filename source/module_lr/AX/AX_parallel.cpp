#ifdef __MPI
#include "AX.h"
#include "module_base/scalapack_connector.h"
#include "module_base/tool_title.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_print.h"
namespace LR
{
    //output: col first, consistent with blas
    // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
    // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
    void cal_AX_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<double>& c,
        const Parallel_2D& pc,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pX,
        double* AX_istate,
        const bool add_on)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_pblas");
        assert(pmat.comm() == pc.comm() && pmat.comm() == pX.comm());
        assert(pmat.blacs_ctxt == pc.blacs_ctxt && pmat.blacs_ctxt == pX.blacs_ctxt);
        assert(pX.get_local_size() > 0);

        const int nks = V_istate.size();

        Parallel_2D pVc;        // for intermediate Vc
        LR_Util::setup_2d_division(pVc, pmat.get_block_size(), naos, nocc, pmat.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            const int ax_start = isk * pX.get_local_size();
            c.fix_k(isk);

            //Vc
            container::Tensor Vc(DAT::DT_DOUBLE, DEV::CpuDevice, { pVc.get_col_size(), pVc.get_row_size() });//row is "inside"(memory contiguity) for pblas
            Vc.zero();

            int i1 = 1;
            int ivirt = nocc + 1;

            char transa = 'N';
            char transb = 'N';
            const double alpha = 1.0;
            const double beta = add_on ? 1.0 : 0.0;
            pdgemm_(&transa, &transb, &naos, &nocc, &naos,
                &alpha, V_istate[isk].data<double>(), &i1, &i1, pmat.desc,
                c.get_pointer(), &i1, &i1, pc.desc,
                &beta, Vc.data<double>(), &i1, &i1, pVc.desc);

            transa = 'T';
            // AX_istate = c ^ TVc
            // descC puts M(nvirt) to row
            pdgemm_(&transa, &transb, &nvirt, &nocc, &naos,
                &alpha, c.get_pointer(), &i1, &ivirt, pc.desc,
                Vc.data<double>(), &i1, &i1, pVc.desc,
                &beta, AX_istate + ax_start, &i1, &i1, pX.desc);

        }
    }

    void cal_AX_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<std::complex<double>>& c,
        const Parallel_2D& pc,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pX,
        std::complex<double>* AX_istate,
        const bool add_on)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_plas");
        assert(pmat.comm() == pc.comm() && pmat.comm() == pX.comm());
        assert(pmat.blacs_ctxt == pc.blacs_ctxt && pmat.blacs_ctxt == pX.blacs_ctxt);
        assert(pX.get_local_size() > 0);

        const int nks = V_istate.size();

        Parallel_2D pVc;        // for intermediate Vc
        LR_Util::setup_2d_division(pVc, pmat.get_block_size(), naos, nocc, pmat.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            const int ax_start = isk * pX.get_local_size();
            c.fix_k(isk);

            //Vc
            container::Tensor Vc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pVc.get_col_size(), pVc.get_row_size() });
            Vc.zero();

            int i1 = 1;
            int ivirt = nocc + 1;

            char transa = 'N';
            char transb = 'N';
            const std::complex<double> alpha(1.0, 0.0);
            const std::complex<double> beta = add_on ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
            pzgemm_(&transa, &transb, &naos, &nocc, &naos,
                &alpha, V_istate[isk].data<std::complex<double>>(), &i1, &i1, pmat.desc,
                c.get_pointer(), &i1, &i1, pc.desc,
                &beta, Vc.data<std::complex<double>>(), &i1, &i1, pVc.desc);

            transa = 'C';
            // AX_istate = c ^ TVc
            // descC puts M(nvirt) to row
            pzgemm_(&transa, &transb, &nvirt, &nocc, &naos,
                &alpha, c.get_pointer(), &i1, &ivirt, pc.desc,
                Vc.data<std::complex<double>>(), &i1, &i1, pVc.desc,
                &beta, AX_istate + ax_start, &i1, &i1, pX.desc);
        }
    }
}
#endif
