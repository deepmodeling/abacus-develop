#ifdef __MPI
#include "CVCX.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
namespace LR
{
    template <>
    void CVCX_occ_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<double, base_device::DEVICE_CPU>& c,
        const Parallel_2D& pc,
        const double* const X_istate,
        const Parallel_2D& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* const AX_istate,
        const bool add_on,
        const double factor)
    {
        ModuleBase::TITLE("hamilt_lrtd", "CVCX_occ_pblas");
        assert(pmat.comm() == pc.comm());
        assert(pmat.comm() == px.comm());
        assert(pmat.blacs_ctxt == pc.blacs_ctxt);
        assert(pmat.blacs_ctxt == px.blacs_ctxt);
        assert(px.get_local_size() > 0);

        const int nks = c.get_nk();
        assert(V_istate.size() == nks);

        Parallel_2D pcv;
        LR_Util::setup_2d_division(pcv, pmat.get_block_size(), nocc, naos, pmat.blacs_ctxt);
        Parallel_2D pcx;
        LR_Util::setup_2d_division(pcx, pmat.get_block_size(), naos, nvirt, pmat.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int start = isk * px.get_local_size();

            const int i1 = 1;
            const int ivirt = nocc + 1;
            const char trans = 'T';
            const char notrans = 'N';  //c is col major
            const double one = 1.0;
            const double zero = 0.0;

            // c^TV[nocc*naos]
            container::Tensor cv(DAT::DT_DOUBLE, DEV::CpuDevice, { pcv.get_col_size(), pcv.get_row_size() });
            pdgemm_(&trans, &notrans, &nocc, &naos, &naos,
                &one, c.get_pointer(), &i1, &i1, pc.desc,
                V_istate[isk].data<double>(), &i1, &i1, pmat.desc,
                &zero, cv.data<double>(), &i1, &i1, pcv.desc);

            // cX^T[naos*nvirt]
            container::Tensor cx(DAT::DT_DOUBLE, DEV::CpuDevice, { pcx.get_col_size(), pcx.get_row_size() });
            pdgemm_(&notrans, &trans, &naos, &nvirt, &nocc,
                &one, c.get_pointer(), &i1, &i1, pc.desc,
                X_istate + start, &i1, &i1, px.desc,
                &zero, cx.data<double>(), &i1, &i1, pcx.desc);

            //AX_istate=[cX^T]^T[c^TV]^T (nvirt major)
            pdgemm_(&trans, &trans, &nvirt, &nocc, &naos,
                &one, cx.data<double>(), &i1, &i1, pcx.desc,
                cv.data<double>(), &i1, &i1, pcv.desc,
                add_on ? &factor : &zero, AX_istate + start, &i1, &i1, px.desc);
        }
    }

    template <>
    void CVCX_occ_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>& c,
        const Parallel_2D& pc,
        const std::complex<double>* const X_istate,
        const Parallel_2D& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const AX_istate,
        const bool add_on,
        const std::complex<double> factor)
    {
        ModuleBase::TITLE("hamilt_lrtd", "CVCX_occ_pblas");
        assert(pmat.comm() == pc.comm());
        assert(pmat.comm() == px.comm());
        assert(pmat.blacs_ctxt == pc.blacs_ctxt);
        assert(pmat.blacs_ctxt == px.blacs_ctxt);
        assert(px.get_local_size() > 0);

        int nks = c.get_nk();
        assert(V_istate.size() == nks);

        Parallel_2D pcv;
        LR_Util::setup_2d_division(pcv, pmat.get_block_size(), nocc, naos, pmat.blacs_ctxt);
        Parallel_2D pcx;
        LR_Util::setup_2d_division(pcx, pmat.get_block_size(), naos, nvirt, pmat.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int start = isk * px.get_local_size();

            const int i1 = 1;
            const int ivirt = nocc + 1;
            const char trans = 'T';
            const char dagger = 'C';
            const char notrans = 'N';  //c is col major
            const std::complex<double> one(1.0, 0.0);
            const std::complex<double> zero(0.0, 0.0);

            // c^TV[nocc*naos]
            container::Tensor cv(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pcv.get_col_size(), pcv.get_row_size() });
            pzgemm_(&dagger, &notrans, &nocc, &naos, &naos,
                &one, c.get_pointer(), &i1, &i1, pc.desc,
                V_istate[isk].data<std::complex<double>>(), &i1, &i1, pmat.desc,
                &zero, cv.data<std::complex<double>>(), &i1, &i1, pcv.desc);

            // cX^T[naos*nvirt]
            container::Tensor cx(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pcx.get_col_size(), pcx.get_row_size() });
            pzgemm_(&notrans, &dagger, &naos, &nvirt, &nocc,
                &one, c.get_pointer(), &i1, &i1, pc.desc,
                X_istate + start, &i1, &i1, px.desc,
                &zero, cx.data<std::complex<double>>(), &i1, &i1, pcx.desc);

            //AX_istate=[cX^T]^T[c^TV]^T (nvirt major)
            pzgemm_(&trans, &trans, &nvirt, &nocc, &naos,
                &one, cx.data<std::complex<double>>(), &i1, &i1, pcx.desc,
                cv.data<std::complex<double>>(), &i1, &i1, pcv.desc,
                add_on ? &factor : &zero, AX_istate + start, &i1, &i1, px.desc);
        }
    }

    template <>
    void CVCX_virt_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<double, base_device::DEVICE_CPU>& c,
        const Parallel_2D& pc,
        const double* const X_istate,
        const Parallel_2D& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* const AX_istate,
        const bool add_on,
        const double factor)
    {
        ModuleBase::TITLE("hamilt_lrtd", "CVCX_virt_pblas");
        assert(pmat.comm() == pc.comm());
        assert(pmat.comm() == px.comm());
        assert(pmat.blacs_ctxt == pc.blacs_ctxt);
        assert(pmat.blacs_ctxt == px.blacs_ctxt);
        assert(px.get_local_size() > 0);

        const int nks = c.get_nk();
        assert(V_istate.size() == nks);

        Parallel_2D pcv;
        LR_Util::setup_2d_division(pcv, pmat.get_block_size(), naos, nvirt, pmat.blacs_ctxt);
        Parallel_2D pcx;
        LR_Util::setup_2d_division(pcx, pmat.get_block_size(), nocc, naos, pmat.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int start = isk * px.get_local_size();

            const int i1 = 1;
            const int ivirt = nocc + 1;
            const char trans = 'T';
            const char dagger = 'C';
            const char notrans = 'N';  //c is col major
            const double one = 1.0;
            const double zero = 0.0;

            // VC[naos*nvirt]
            container::Tensor cv(DAT::DT_DOUBLE, DEV::CpuDevice, { pcv.get_col_size(), pcv.get_row_size() });
            pdgemm_(&notrans, &notrans, &naos, &nvirt, &naos,
                &one, V_istate[isk].data<double>(), &i1, &i1, pmat.desc,
                c.get_pointer(), &i1, &ivirt, pc.desc,
                &zero, cv.data<double>(), &i1, &i1, pcv.desc);

            // X^TC^T[nocc*naos]
            container::Tensor cx(DAT::DT_DOUBLE, DEV::CpuDevice, { pcx.get_col_size(), pcx.get_row_size() });
            pdgemm_(&dagger, &dagger, &nocc, &naos, &nvirt,
                &one, X_istate + start, &i1, &i1, px.desc,
                c.get_pointer(), &i1, &ivirt, pc.desc,
                &zero, cx.data<double>(), &i1, &i1, pcx.desc);

            //AX_istate=[VC]^T[X^TC^T]^T (nvirt major)
            pdgemm_(&trans, &trans, &nvirt, &nocc, &naos,
                &one, cv.data<double>(), &i1, &i1, pcv.desc,
                cx.data<double>(), &i1, &i1, pcx.desc,
                add_on ? &factor : &zero, AX_istate + start, &i1, &i1, px.desc);
        }
    }

    template <>
    void CVCX_virt_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>& c,
        const Parallel_2D& pc,
        const std::complex<double>* const X_istate,
        const Parallel_2D& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const AX_istate,
        const bool add_on,
        const std::complex<double> factor)
    {
        ModuleBase::TITLE("hamilt_lrtd", "CVCX_virt_pblas");
        assert(pmat.comm() == pc.comm());
        assert(pmat.comm() == px.comm());
        assert(pmat.blacs_ctxt == pc.blacs_ctxt);
        assert(pmat.blacs_ctxt == px.blacs_ctxt);
        assert(px.get_local_size() > 0);

        const int nks = c.get_nk();
        assert(V_istate.size() == nks);

        Parallel_2D pcv;
        LR_Util::setup_2d_division(pcv, pmat.get_block_size(), naos, nvirt, pmat.blacs_ctxt);
        Parallel_2D pcx;
        LR_Util::setup_2d_division(pcx, pmat.get_block_size(), nocc, naos, pmat.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int start = isk * px.get_local_size();

            const int i1 = 1;
            const int ivirt = nocc + 1;
            const char trans = 'T';
            const char dagger = 'C';
            const char notrans = 'N';  //c is col major
            const std::complex<double> one(1.0, 0.0);
            const std::complex<double> zero(0.0, 0.0);

            // VC[naos*nvirt]
            container::Tensor cv(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pcv.get_col_size(), pcv.get_row_size() });
            pzgemm_(&notrans, &notrans, &naos, &nvirt, &naos,
                &one, V_istate[isk].data<std::complex<double>>(), &i1, &i1, pmat.desc,
                c.get_pointer(), &i1, &ivirt, pc.desc,
                &zero, cv.data<std::complex<double>>(), &i1, &i1, pcv.desc);

            // X^TC^T[nocc*naos]
            container::Tensor cx(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pcx.get_col_size(), pcx.get_row_size() });
            pzgemm_(&dagger, &dagger, &nocc, &naos, &nvirt,
                &one, X_istate + start, &i1, &i1, px.desc,
                c.get_pointer(), &i1, &ivirt, pc.desc,
                &zero, cx.data<std::complex<double>>(), &i1, &i1, pcx.desc);

            //AX_istate=[VC]^T[X^TC^T]^T (nvirt major)
            pzgemm_(&trans, &trans, &nvirt, &nocc, &naos,
                &one, cv.data<std::complex<double>>(), &i1, &i1, pcv.desc,
                cx.data<std::complex<double>>(), &i1, &i1, pcx.desc,
                add_on ? &factor : &zero, AX_istate + start, &i1, &i1, px.desc);
        }
    }
}
#endif