#include "cal_ldos.h"

#include "cube_io.h"
#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"

#include <type_traits>

namespace ModuleIO
{
template <typename T>
void Cal_ldos<T>::cal_ldos_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                              const psi::Psi<std::complex<double>>& psi,
                              const Parallel_Grid& pgrid,
                              const UnitCell& ucell)
{
    // energy range for ldos (efermi as reference)
    const double emin = PARAM.inp.stm_bias < 0 ? PARAM.inp.stm_bias : 0;
    const double emax = PARAM.inp.stm_bias > 0 ? PARAM.inp.stm_bias : 0;

    std::vector<double> ldos(pelec->charge->nrxx);
    std::vector<std::complex<double>> wfcr(pelec->basis->nrxx);

    for (int ik = 0; ik < pelec->klist->get_nks(); ++ik)
    {
        psi.fix_k(ik);
        const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);
        int nbands = psi.get_nbands();

        for (int ib = 0; ib < nbands; ib++)
        {
            pelec->basis->recip2real(&psi(ib, 0), wfcr.data(), ik);
            const double eigenval = (pelec->ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
            if (eigenval >= emin && eigenval <= emax)
            {
                for (int ir = 0; ir < pelec->basis->nrxx; ir++)
                    ldos[ir] += pelec->klist->wk[ik] * norm(wfcr[ir]);
            }
        }
    }

    std::stringstream fn;
    fn << PARAM.globalv.global_out_dir << "LDOS_" << PARAM.inp.stm_bias << "eV"
       << ".cube";

    const int precision = PARAM.inp.out_ldos[1];
    ModuleIO::write_vdata_palgrid(pgrid, ldos.data(), 0, PARAM.inp.nspin, 0, fn.str(), 0, &ucell, precision, 0);
}

#ifdef __LCAO
template <typename T>
void Cal_ldos<T>::cal_ldos_lcao(const elecstate::ElecStateLCAO<T>* pelec,
                                const psi::Psi<T>& psi,
                                const Parallel_Grid& pgrid,
                                const UnitCell& ucell)
{
    // energy range for ldos (efermi as reference)
    const double emin = PARAM.inp.stm_bias < 0 ? PARAM.inp.stm_bias : 0;
    const double emax = PARAM.inp.stm_bias > 0 ? PARAM.inp.stm_bias : 0;

    // calulate dm-like
    const int nbands_local = psi.get_nbands();
    const int nbasis_local = psi.get_nbasis();

    // psi.T * wk * psi.conj()
    // result[ik](iw1,iw2) = \sum_{ib} psi[ik](ib,iw1).T * wk(k) * psi[ik](ib,iw2).conj()
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);

        // T* dmk_pointer = DM.get_DMK_pointer(ik);

        psi::Psi<T> wk_psi(1, psi.get_nbands(), psi.get_nbasis(), psi.get_nbasis(), true);
        const T* ppsi = psi.get_pointer();
        T* pwk_psi = wk_psi.get_pointer();

        // #ifdef _OPENMP
        // #pragma omp parallel for schedule(static, 1024)
        // #endif
        //         for (int i = 0; i < wk_psi.size(); ++i)
        //         {
        //             pwk_psi[i] = my_conj(ppsi[i]);
        //         }

        //         int ib_global = 0;
        //         for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        //         {
        //             while (ib_local != ParaV->global2local_col(ib_global))
        //             {
        //                 ++ib_global;
        //                 if (ib_global >= wg.nc)
        //                 {
        //                     ModuleBase::WARNING_QUIT("cal_ldos", "please check global2local_col!");
        //                 }
        //             }

        //             const double eigenval = (pelec->ekb(ik, ib_global) - efermi) * ModuleBase::Ry_to_eV;
        //             if (eigenval >= emin && eigenval <= emax)
        //             {
        //                 for (int ir = 0; ir < pelec->basis->nrxx; ir++)
        //                     ldos[ir] += pelec->klist->wk[ik] * norm(wfcr[ir]);
        //             }

        //             double* wg_wfc_pointer = &(wk_psi(0, ib_local, 0));
        //             BlasConnector::scal(nbasis_local, pelec->klist->wk[ik], wg_wfc_pointer, 1);
        //         }

        //         // C++: dm(iw1,iw2) = psi(ib,iw1).T * wk_psi(ib,iw2)
        // #ifdef __MPI
        //         psiMulPsiMpi(wk_psi, psi, dmk_pointer, ParaV->desc_wfc, ParaV->desc);
        // #else
        //         psiMulPsi(wk_psi, psi, dmk_pointer);
        // #endif
    }
}

double my_conj(double x)
{
    return x;
}

std::complex<double> my_conj(const std::complex<double>& z)
{
    return {z.real(), -z.imag()};
}

#endif

template class Cal_ldos<double>;               // Gamma_only case
template class Cal_ldos<std::complex<double>>; // multi-k case
} // namespace elecstate
