#include "force_rt_overlap.h"

#include "source_base/basis_flattened_indexing.h"
#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_reduce.h"
#include "td_folding.h"
#include "td_info.h"
template <>
void cal_foverlap_rt(ModuleBase::matrix& foverlap,
                     const module_dm::Setup_DM<std::complex<double>>& dmat,
                     hamilt::Hamilt<std::complex<double>>* p_hamilt,
                     const K_Vectors& kv,
                     const BasisFlattenedIndexing& basis,
                     const BlockCyclicDistributionInformation& layout,
                     UnitCell& ucell,
                     TD_info& td_state)
{
#ifdef __MPI
    const int nlocal = layout.nrglb();
    assert(nlocal >= 0);

    TD_info* td_info = &td_state;
    // get dS/dR_{x,y,z}
    std::vector<hamilt::HContainer<double>*> dsxr = td_info->get_grad_overlap();
    // allocate matrix
    const std::size_t nloc = layout.size();
    const int nrow = layout.nrloc();
    std::vector<std::complex<double>> Htmp(nloc);
    std::vector<std::complex<double>> Sinv(nloc);
    std::vector<std::complex<double>> dsxk(nloc);
    std::vector<std::complex<double>> pdsxk(nloc);
    std::vector<std::complex<double>> tmp1(nloc);
    std::vector<std::complex<double>> tmp2(nloc);
    std::vector<std::complex<double>> tmp3(nloc);
    std::vector<std::vector<std::complex<double>>> tmp_out(3, std::vector<std::complex<double>>(nloc, 0.0));
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        p_hamilt->updateHk(ik);
        // get dmk
        std::complex<double>* tmp_dmk = dmat.dm->get_dmk_ptr(ik);

        ModuleBase::GlobalFunc::ZEROS(Htmp.data(), nloc);
        ModuleBase::GlobalFunc::ZEROS(Sinv.data(), nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp1.data(), nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp2.data(), nloc);

        const int inc = 1;

        hamilt::MatrixBlock<std::complex<double>> h_mat;
        hamilt::MatrixBlock<std::complex<double>> s_mat;
        // get Hk Sk
        p_hamilt->matrix(h_mat, s_mat);
        BlasConnector::copy(nloc, h_mat.p, inc, Htmp.data(), inc);
        BlasConnector::copy(nloc, s_mat.p, inc, Sinv.data(), inc);

        std::vector<int> ipiv(nloc, 0);
        int info = 0;
        const int one_int = 1;
        const std::complex<double> one_complex = {1.0, 0.0};
        const std::complex<double> two_complex = {2.0, 0.0};
        const std::complex<double> mone_complex = {-1.0, 0.0};
        const std::complex<double> zero_complex = {0.0, 0.0};

        ScalapackConnector::getrf(nlocal, nlocal, Sinv.data(), one_int, one_int, layout.desc(), ipiv.data(), &info);

        int lwork = -1;
        int liwork = -1;

        // if lwork == -1, then the size of work is (at least) of length 1.
        std::vector<std::complex<double>> work(1, 0);

        // if liwork = -1, then the size of iwork is (at least) of length 1.
        std::vector<int> iwork(1, 0);

        ScalapackConnector::getri(nlocal,
                                  Sinv.data(),
                                  one_int,
                                  one_int,
                                  layout.desc(),
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
                                  Sinv.data(),
                                  one_int,
                                  one_int,
                                  layout.desc(),
                                  ipiv.data(),
                                  work.data(),
                                  &lwork,
                                  iwork.data(),
                                  &liwork,
                                  &info);

        const char N_char = 'N';
        const char T_char = 'T';
        const char C_char = 'C';

        ScalapackConnector::gemm(T_char,
                                 C_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp_dmk,
                                 one_int,
                                 one_int,
                                 layout.desc(),
                                 Htmp.data(),
                                 one_int,
                                 one_int,
                                 layout.desc(),
                                 zero_complex,
                                 tmp1.data(),
                                 one_int,
                                 one_int,
                                 layout.desc());

        ScalapackConnector::gemm(N_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp1.data(),
                                 one_int,
                                 one_int,
                                 layout.desc(),
                                 Sinv.data(),
                                 one_int,
                                 one_int,
                                 layout.desc(),
                                 zero_complex,
                                 tmp2.data(),
                                 one_int,
                                 one_int,
                                 layout.desc());
        for (int dir = 0; dir < 3; dir++)
        {
            ModuleBase::GlobalFunc::ZEROS(dsxk.data(), nloc);
            ModuleBase::GlobalFunc::ZEROS(pdsxk.data(), nloc);
            ModuleBase::GlobalFunc::ZEROS(tmp3.data(), nloc);
            module_rt::folding_HR_td(*dsxr[dir], dsxk.data(), kv.kvec_d[ik], td_info->endpoint_A(), td_info->get_phase_hybrid(), nrow, 1);
            const ModuleBase::Vector3<double> electric_field_ry = td_info->field().efield_ha() * 2.0;
            module_rt::folding_partial_dot(*dsxr[dir],
                                           pdsxk.data(),
                                           kv.kvec_d[ik],
                                           nrow,
                                           1,
                                           &ucell,
                                           td_info->get_phase_hybrid(),
                                           td_info->endpoint_A(),
                                           electric_field_ry);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     two_complex,
                                     tmp2.data(),
                                     one_int,
                                     one_int,
                                     layout.desc(),
                                     dsxk.data(),
                                     one_int,
                                     one_int,
                                     layout.desc(),
                                     one_complex,
                                     tmp_out[dir].data(),
                                     one_int,
                                     one_int,
                                     layout.desc());

            ScalapackConnector::gemm(T_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     mone_complex,
                                     tmp_dmk,
                                     one_int,
                                     one_int,
                                     layout.desc(),
                                     pdsxk.data(),
                                     one_int,
                                     one_int,
                                     layout.desc(),
                                     one_complex,
                                     tmp3.data(),
                                     one_int,
                                     one_int,
                                     layout.desc());

            ScalapackConnector::geadd(N_char,
                                      nlocal,
                                      nlocal,
                                      one_complex,
                                      tmp3.data(),
                                      one_int,
                                      one_int,
                                      layout.desc(),
                                      one_complex,
                                      tmp_out[dir].data(),
                                      one_int,
                                      one_int,
                                      layout.desc());
        }
    }
#pragma omp parallel for
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        double* force_tmp1 = &foverlap(iat, 0);
        const BasisRange range = basis.range(iat);
        std::vector<int> row_positions;
        std::vector<int> col_positions;
        for (int iglb = static_cast<int>(range.ibas0); iglb < static_cast<int>(range.ibas1); ++iglb)
        {
            const int irloc = layout.irloc(iglb);
            if (irloc >= 0)
            {
                row_positions.push_back(irloc);
            }
            const int icloc = layout.icloc(iglb);
            if (icloc >= 0)
            {
                col_positions.push_back(icloc);
            }
        }
        for (const int irloc : row_positions)
        {
            for (const int icloc : col_positions)
            {
                if (layout.irglb(irloc) == layout.icglb(icloc))
                {
                    const std::size_t index = static_cast<std::size_t>(irloc)
                                              + static_cast<std::size_t>(icloc) * static_cast<std::size_t>(nrow);
                    for (int dir = 0; dir < 3; dir++)
                    {
                        force_tmp1[dir] += (tmp_out[dir][index]).real();
                    }
                }
            }
        }
    }
    Parallel_Reduce::reduce_all(foverlap.c, foverlap.nr * foverlap.nc);
    return;
#endif
}
template <>
void cal_foverlap_rt(ModuleBase::matrix& foverlap,
                     const module_dm::Setup_DM<double>& dmat,
                     hamilt::Hamilt<double>* p_hamilt,
                     const K_Vectors& kv,
                     const BasisFlattenedIndexing& basis,
                     const BlockCyclicDistributionInformation& layout,
                     UnitCell& ucell,
                     TD_info& td_state)
{
    return;
}
