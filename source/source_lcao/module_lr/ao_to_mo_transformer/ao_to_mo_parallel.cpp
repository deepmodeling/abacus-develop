#ifdef __MPI
#include "ao_to_mo.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
namespace LR
{
    //output: col first, consistent with blas
    // coeff: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
    // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
    template<>
    void ao_to_mo_pblas(
        const std::vector<container::Tensor>& mat_ao,
        const Parallel_2D& pmat_ao,
        const psi::Psi<double>& coeff,
        const Parallel_2D& pcoeff,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pmat_mo,
        double* mat_mo,
        const bool add_on,
        const LR_Util::MO_TYPE type)
    {
        ModuleBase::TITLE("LR", "ao_to_mo_pblas");
        assert(pmat_ao.comm() == pcoeff.comm() && pmat_ao.comm() == pmat_mo.comm());
        assert(pmat_ao.blacs_ctxt == pcoeff.blacs_ctxt && pmat_ao.blacs_ctxt == pmat_mo.blacs_ctxt);
        assert(pmat_mo.get_local_size() > 0);

        const int nks = mat_ao.size();
        int nmo1_set, nmo2_set, imo1_set, imo2_set;
        LR_Util::set_dim(type, nocc, nvirt, nmo1_set, nmo2_set, imo1_set, imo2_set);
        const int nmo1 = nmo1_set;
        const int nmo2 = nmo2_set;
        const int imo1 = imo1_set + 1;
        const int imo2 = imo2_set + 1;
        const int i1 = 1;

        Parallel_2D pVc;        // for intermediate Vc
        LR_Util::setup_2d_division(pVc, pmat_ao.get_block_size(), naos, nmo1, pmat_ao.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            const int start = isk * pmat_mo.get_local_size();
            coeff.fix_k(isk);

            //Vc
            container::Tensor Vc(DAT::DT_DOUBLE, DEV::CpuDevice, { pVc.get_col_size(), pVc.get_row_size() });//row is "inside"(memory contiguity) for pblas
            Vc.zero();

            char transa = 'N';
            char transb = 'N';
            const double alpha = 1.0;
            const double beta = add_on ? 1.0 : 0.0;
            ScalapackConnector::gemm(transa, transb, naos, nmo1, naos,
                alpha, mat_ao[isk].data<double>(), i1, i1, pmat_ao.desc,
                coeff.get_pointer(), i1, imo1, pcoeff.desc,
                beta, Vc.data<double>(), i1, i1, pVc.desc);

            transa = 'T';
            // mat_mo = c ^ TVc
            // descC puts M(nvirt) to row
            ScalapackConnector::gemm(transa, transb, nmo2, nmo1, naos,
                alpha, coeff.get_pointer(), i1, imo2, pcoeff.desc,
                Vc.data<double>(), i1, i1, pVc.desc,
                beta, mat_mo + start, i1, i1, pmat_mo.desc);

        }
    }

    template<>
    void ao_to_mo_pblas(
        const std::vector<container::Tensor>& mat_ao,
        const Parallel_2D& pmat_ao,
        const psi::Psi<std::complex<double>>& coeff,
        const Parallel_2D& pcoeff,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pmat_mo,
        std::complex<double>* const mat_mo,
        const bool add_on,
        const LR_Util::MO_TYPE type)
    {
        ModuleBase::TITLE("LR", "ao_to_mo_pblas");
        assert(pmat_ao.comm() == pcoeff.comm() && pmat_ao.comm() == pmat_mo.comm());
        assert(pmat_ao.blacs_ctxt == pcoeff.blacs_ctxt && pmat_ao.blacs_ctxt == pmat_mo.blacs_ctxt);
        assert(pmat_mo.get_local_size() > 0);

        const int nks = mat_ao.size();
        int nmo1_set, nmo2_set, imo1_set, imo2_set;
        LR_Util::set_dim(type, nocc, nvirt, nmo1_set, nmo2_set, imo1_set, imo2_set);
        const int nmo1 = nmo1_set;
        const int nmo2 = nmo2_set;
        const int imo1 = imo1_set + 1;
        const int imo2 = imo2_set + 1;
        const int i1 = 1;

        Parallel_2D pVc;        // for intermediate Vc
        LR_Util::setup_2d_division(pVc, pmat_ao.get_block_size(), naos, nmo1, pmat_ao.blacs_ctxt);
        for (int isk = 0;isk < nks;++isk)
        {
            const int start = isk * pmat_mo.get_local_size();
            coeff.fix_k(isk);

            //Vc
            container::Tensor Vc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { pVc.get_col_size(), pVc.get_row_size() });
            Vc.zero();

            char transa = 'N';
            char transb = 'N';
            const std::complex<double> alpha(1.0, 0.0);
            const std::complex<double> beta = add_on ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
            ScalapackConnector::gemm(transa, transb, naos, nmo1, naos,
                alpha, mat_ao[isk].data<std::complex<double>>(), i1, i1, pmat_ao.desc,
                coeff.get_pointer(), i1, imo1, pcoeff.desc,
                beta, Vc.data<std::complex<double>>(), i1, i1, pVc.desc);

            transa = 'C';
            // mat_mo = c ^ TVc
            // descC puts M(nvirt) to row
            ScalapackConnector::gemm(transa, transb, nmo2, nmo1, naos,
                alpha, coeff.get_pointer(), i1, imo2, pcoeff.desc,
                Vc.data<std::complex<double>>(), i1, i1, pVc.desc,
                beta, mat_mo + start, i1, i1, pmat_mo.desc);
        }
    }
}
#endif
