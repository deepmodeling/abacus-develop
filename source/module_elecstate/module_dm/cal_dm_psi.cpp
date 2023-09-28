#include "cal_dm_psi.h"

#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_psi/psi.h"

namespace elecstate
{

// for Gamma-Only case where DMK is double
void cal_dm_psi(const Parallel_Orbitals* ParaV,
                       const ModuleBase::matrix& wg,
                       const psi::Psi<double>& wfc,
                       elecstate::DensityMatrix<double, double>& DM)
{
    ModuleBase::TITLE("elecstate", "cal_dm");
    ModuleBase::timer::tick("elecstate", "cal_dm");

    // dm.resize(wfc.get_nk(), ParaV->ncol, ParaV->nrow);
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // dm = wfc.T * wg * wfc.conj()
    // dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        double* dmk_pointer = DM.get_DMK_pointer(ik);
        wfc.fix_k(ik);
        // dm.fix_k(ik);
        // dm[ik].create(ParaV->ncol, ParaV->nrow);
        //  wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi::Psi<double> wg_wfc(wfc, 1);

        int ib_global = 0;
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            while (ib_local != ParaV->global2local_col(ib_global))
            {
                ++ib_global;
                if (ib_global >= wg.nc)
                {
                    break;
                    ModuleBase::WARNING_QUIT("ElecStateLCAO::cal_dm", "please check global2local_col!");
                }
            }
            if (ib_global >= wg.nc)
                continue;
            const double wg_local = wg(ik, ib_global);
            double* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dmk_pointer, ParaV->desc_wfc, ParaV->desc);
#else
        psiMulPsi(wg_wfc, wfc, dmk_pointer);
#endif
    }
    ModuleBase::timer::tick("elecstate", "cal_dm");

    return;
}

void cal_dm_psi(const Parallel_Orbitals* ParaV,
                       const ModuleBase::matrix& wg,
                       const psi::Psi<std::complex<double>>& wfc,
                       elecstate::DensityMatrix<std::complex<double>, double>& DM)
{
    ModuleBase::TITLE("elecstate", "cal_dm");
    ModuleBase::timer::tick("elecstate", "cal_dm");

    // dm.resize(wfc.get_nk(), ParaV->ncol, ParaV->nrow);
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // dm = wfc.T * wg * wfc.conj()
    // dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        std::complex<double>* dmk_pointer = DM.get_DMK_pointer(ik);
        // dm.fix_k(ik);
        //dm[ik].create(ParaV->ncol, ParaV->nrow);
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi::Psi<std::complex<double>> wg_wfc(1, wfc.get_nbands(), wfc.get_nbasis(), nullptr);
        const std::complex<double>* pwfc = wfc.get_pointer();
        std::complex<double>* pwg_wfc = wg_wfc.get_pointer();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int i = 0; i < wg_wfc.size(); ++i)
        {
            pwg_wfc[i] = conj(pwfc[i]);
        }

        int ib_global = 0;
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            while (ib_local != ParaV->global2local_col(ib_global))
            {
                ++ib_global;
                if (ib_global >= wg.nc)
                {
                    break;
                    ModuleBase::WARNING_QUIT("ElecStateLCAO::cal_dm", "please check global2local_col!");
                }
            }
            if (ib_global >= wg.nc)
                continue;
            const double wg_local = wg(ik, ib_global);
            std::complex<double>* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dmk_pointer, ParaV->desc_wfc, ParaV->desc);
#else
        psiMulPsi(wg_wfc, wfc, dmk_pointer);
#endif
    }

    ModuleBase::timer::tick("elecstate", "cal_dm");
    return;
}

#ifdef __MPI
void psiMulPsiMpi(const psi::Psi<double>& psi1,
                         const psi::Psi<double>& psi2,
                         double* dm_out,
                         const int* desc_psi,
                         const int* desc_dm)
{
    ModuleBase::timer::tick("psiMulPsiMpi", "pdgemm");
    const double one_float = 1.0, zero_float = 0.0;
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    pdgemm_(&N_char,
            &T_char,
            &nlocal,
            &nlocal,
            &nbands,
            &one_float,
            psi1.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            psi2.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            &zero_float,
            dm_out,
            &one_int,
            &one_int,
            desc_dm);
    ModuleBase::timer::tick("psiMulPsiMpi", "pdgemm");
}

void psiMulPsiMpi(const psi::Psi<std::complex<double>>& psi1,
                         const psi::Psi<std::complex<double>>& psi2,
                         std::complex<double>* dm_out,
                         const int* desc_psi,
                         const int* desc_dm)
{
    ModuleBase::timer::tick("psiMulPsiMpi", "pdgemm");
    const std::complex<double> one_complex = {1.0, 0.0}, zero_complex = {0.0, 0.0};
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    pzgemm_(&N_char,
            &T_char,
            &nlocal,
            &nlocal,
            &nbands,
            &one_complex,
            psi1.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            psi2.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            &zero_complex,
            dm_out,
            &one_int,
            &one_int,
            desc_dm);
    ModuleBase::timer::tick("psiMulPsiMpi", "pdgemm");
}

#else
void psiMulPsi(const psi::Psi<double>& psi1, const psi::Psi<double>& psi2, double* dm_out)
{
    const double one_float = 1.0, zero_float = 0.0;
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    dgemm_(&N_char,
           &T_char,
           &nlocal,
           &nlocal,
           &nbands,
           &one_float,
           psi1.get_pointer(),
           &nlocal,
           psi2.get_pointer(),
           &nlocal,
           &zero_float,
           dm_out,
           &nlocal);
}

void psiMulPsi(const psi::Psi<std::complex<double>>& psi1,
                      const psi::Psi<std::complex<double>>& psi2,
                      std::complex<double>* dm_out)
{
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    const std::complex<double> one_complex = {1.0, 0.0}, zero_complex = {0.0, 0.0};
    zgemm_(&N_char,
           &T_char,
           &nlocal,
           &nlocal,
           &nbands,
           &one_complex,
           psi1.get_pointer(),
           &nlocal,
           psi2.get_pointer(),
           &nlocal,
           &zero_complex,
           dm_out,
           &nlocal);
}
#endif

} // namespace elecstate
