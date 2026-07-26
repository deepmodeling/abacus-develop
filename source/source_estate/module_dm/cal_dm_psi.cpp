#include "cal_dm_psi.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_psi/psi.h"

namespace elecstate
{

static void cast_down_to_float(const std::complex<double>* src, std::complex<float>* dst, std::size_t n)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (std::size_t i = 0; i < n; ++i)
    {
        dst[i] = static_cast<std::complex<float>>(src[i]);
    }
}

static void cast_up_to_double(const std::complex<float>* src, std::complex<double>* dst, std::size_t n)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (std::size_t i = 0; i < n; ++i)
    {
        dst[i] = static_cast<std::complex<double>>(src[i]);
    }
}

void cal_dm_psi(const Parallel_Orbitals* ParaV,
                       const ModuleBase::matrix& wg,
                       const psi::Psi<double>& wfc,
                       elecstate::DensityMatrix<double, double>& DM)
{
    ModuleBase::TITLE("elecstate", "cal_dm_psi");
    ModuleBase::timer::start("elecstate", "cal_dm_psi");

    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        double* dmk_pointer = DM.get_DMK_pointer(ik);
        wfc.fix_k(ik);

        psi::Psi<double> wg_wfc(1, wfc.get_nbands(), 
          wfc.get_nbasis(), wfc.get_nbasis(), true);

        wg_wfc.set_all_psi(wfc.get_pointer(), wg_wfc.size());

        int ib_global = 0;
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            while (ib_local != ParaV->global2local_col(ib_global))
            {
                ++ib_global;
                if (ib_global >= wg.nc)
                {
                    break;
                }
            }
			if (ib_global >= wg.nc)
			{
				continue;
			}
            const double wg_local = wg(ik, ib_global);

            double* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dmk_pointer, ParaV->desc_wfc, ParaV->desc);
#else
        psiMulPsi(wg_wfc, wfc, dmk_pointer);
#endif
    }
    ModuleBase::timer::end("elecstate", "cal_dm_psi");

    return;
}
template <typename TR>
void cal_dm_psi(const Parallel_Orbitals* ParaV,
                       const ModuleBase::matrix& wg,
                       const psi::Psi<std::complex<double>>& wfc,
                       elecstate::DensityMatrix<std::complex<double>, TR>& DM)
{
    ModuleBase::TITLE("elecstate", "cal_dm_psi");
    ModuleBase::timer::start("elecstate", "cal_dm_psi");

    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        std::complex<double>* dmk_pointer = DM.get_DMK_pointer(ik);
        psi::Psi<std::complex<double>> wg_wfc(1, 
                                              wfc.get_nbands(), 
                                              wfc.get_nbasis(),
                                              wfc.get_nbasis(),
                                              true);
        
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
                }
            }
			if (ib_global >= wg.nc)
			{
				continue;
			}
			const double wg_local = wg(ik, ib_global);
            std::complex<double>* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

#ifdef __MPI
        if (PARAM.inp.ks_solver == "cg_in_lcao")
        {
            psiMulPsi(wg_wfc, wfc, dmk_pointer);
		} 
		else 
		{
            psiMulPsiMpi(wg_wfc, wfc, dmk_pointer, ParaV->desc_wfc, ParaV->desc);
        }
#else
        psiMulPsi(wg_wfc, wfc, dmk_pointer);
#endif
    }

    ModuleBase::timer::end("elecstate", "cal_dm_psi");
    return;
}

void cal_dm_psi_mixed(const Parallel_Orbitals* ParaV,
                             const ModuleBase::matrix& wg,
                             const psi::Psi<std::complex<double>>& wfc,
                             elecstate::DensityMatrix<std::complex<double>, double>& DM)
{
    ModuleBase::TITLE("elecstate", "cal_dm_psi_mixed");
    ModuleBase::timer::start("elecstate", "cal_dm_psi_mixed");

    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        std::complex<double>* dmk_pointer = DM.get_DMK_pointer(ik);
        psi::Psi<std::complex<double>> wg_wfc(1,
                                              wfc.get_nbands(),
                                              wfc.get_nbasis(),
                                              wfc.get_nbasis(),
                                              true);

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
                }
            }
            if (ib_global >= wg.nc)
            {
                continue;
            }
            const double wg_local = wg(ik, ib_global);
            std::complex<double>* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

#ifdef __MPI
        if (PARAM.inp.ks_solver == "cg_in_lcao")
        {
            psiMulPsiMixed(wg_wfc, wfc, dmk_pointer);
        }
        else
        {
            psiMulPsiMpiMixed(wg_wfc, wfc, dmk_pointer, ParaV->desc_wfc, ParaV->desc, ParaV->nloc);
        }
#else
        psiMulPsiMixed(wg_wfc, wfc, dmk_pointer);
#endif
    }

    ModuleBase::timer::end("elecstate", "cal_dm_psi_mixed");
    return;
}

#ifdef __MPI
void psiMulPsiMpi(const psi::Psi<double>& psi1,
                         const psi::Psi<double>& psi2,
                         double* dm_out,
                         const int* desc_psi,
                         const int* desc_dm)
{
    ModuleBase::timer::start("psiMulPsiMpi", "pdgemm");
    const double one_float = 1.0, zero_float = 0.0;
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];

    ScalapackConnector::gemm(N_char,
            T_char,
            nlocal,
            nlocal,
            nbands,
            one_float,
            psi1.get_pointer(),
            one_int,
            one_int,
            desc_psi,
            psi2.get_pointer(),
            one_int,
            one_int,
            desc_psi,
            zero_float,
            dm_out,
            one_int,
            one_int,
            desc_dm);
    ModuleBase::timer::end("psiMulPsiMpi", "pdgemm");
}

void psiMulPsiMpi(const psi::Psi<std::complex<double>>& psi1,
                         const psi::Psi<std::complex<double>>& psi2,
                         std::complex<double>* dm_out,
                         const int* desc_psi,
                         const int* desc_dm)
{
    ModuleBase::timer::start("psiMulPsiMpi", "pzgemm");
    const std::complex<double> one_complex = {1.0, 0.0}, zero_complex = {0.0, 0.0};
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    ScalapackConnector::gemm(N_char,
            T_char,
            nlocal,
            nlocal,
            nbands,
            one_complex,
            psi1.get_pointer(),
            one_int,
            one_int,
            desc_psi,
            psi2.get_pointer(),
            one_int,
            one_int,
            desc_psi,
            zero_complex,
            dm_out,
            one_int,
            one_int,
            desc_dm);
    ModuleBase::timer::end("psiMulPsiMpi", "pzgemm");
}

void psiMulPsiMpiMixed(const psi::Psi<std::complex<double>>& psi1,
                              const psi::Psi<std::complex<double>>& psi2,
                              std::complex<double>* dm_out,
                              const int* desc_psi,
                              const int* desc_dm,
                              int nloc_dm)
{
    ModuleBase::timer::start("psiMulPsiMpiMixed", "pcgemm");
    const std::complex<float> one_f = {1.0f, 0.0f}, zero_f = {0.0f, 0.0f};
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];

    std::size_t nloc_psi = psi1.size();
    std::vector<std::complex<float>> psi1_f(nloc_psi);
    cast_down_to_float(psi1.get_pointer(), psi1_f.data(), nloc_psi);
    std::vector<std::complex<float>> psi2_f(psi2.size());
    cast_down_to_float(psi2.get_pointer(), psi2_f.data(), psi2.size());

    std::vector<std::complex<float>> dm_f(nloc_dm, zero_f);

    ScalapackConnector::gemm(N_char,
            T_char,
            nlocal,
            nlocal,
            nbands,
            one_f,
            psi1_f.data(),
            one_int,
            one_int,
            desc_psi,
            psi2_f.data(),
            one_int,
            one_int,
            desc_psi,
            zero_f,
            dm_f.data(),
            one_int,
            one_int,
            desc_dm);
    cast_up_to_double(dm_f.data(), dm_out, nloc_dm);
    ModuleBase::timer::end("psiMulPsiMpiMixed", "pcgemm");
}

#endif

void psiMulPsi(const psi::Psi<double>& psi1, const psi::Psi<double>& psi2, double* dm_out)
{
    const double one_float = 1.0, zero_float = 0.0;
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    BlasConnector::gemm_cm(N_char,
           T_char,
           nlocal,
           nlocal,
           nbands,
           one_float,
           psi1.get_pointer(),
           nlocal,
           psi2.get_pointer(),
           nlocal,
           zero_float,
           dm_out,
           nlocal);
}

void psiMulPsi(const psi::Psi<std::complex<double>>& psi1,
                      const psi::Psi<std::complex<double>>& psi2,
                      std::complex<double>* dm_out)
{
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    BlasConnector::gemm_cm(N_char,
           T_char,
           nlocal,
           nlocal,
           nbands,
           one_complex,
           psi1.get_pointer(),
           nlocal,
           psi2.get_pointer(),
           nlocal,
           zero_complex,
           dm_out,
           nlocal);
}

void psiMulPsiMixed(const psi::Psi<std::complex<double>>& psi1,
                          const psi::Psi<std::complex<double>>& psi2,
                          std::complex<double>* dm_out)
{
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    const std::complex<float> one_f = {1.0f, 0.0f};
    const std::complex<float> zero_f = {0.0f, 0.0f};

    std::size_t sz = psi1.size();
    std::vector<std::complex<float>> psi1_f(sz);
    cast_down_to_float(psi1.get_pointer(), psi1_f.data(), sz);
    std::vector<std::complex<float>> psi2_f(psi2.size());
    cast_down_to_float(psi2.get_pointer(), psi2_f.data(), psi2.size());

    std::vector<std::complex<float>> dm_f(nlocal * nlocal, zero_f);
    BlasConnector::gemm_cm(N_char,
           T_char,
           nlocal,
           nlocal,
           nbands,
           one_f,
           psi1_f.data(),
           nlocal,
           psi2_f.data(),
           nlocal,
           zero_f,
           dm_f.data(),
           nlocal);
    cast_up_to_double(dm_f.data(), dm_out, static_cast<std::size_t>(nlocal) * nlocal);
}

} // namespace elecstate

// Explicit template instantiations
namespace elecstate
{
template void cal_dm_psi<double>(const Parallel_Orbitals*, const ModuleBase::matrix&, const psi::Psi<std::complex<double>>&, DensityMatrix<std::complex<double>, double>&);
template void cal_dm_psi<std::complex<double>>(const Parallel_Orbitals*, const ModuleBase::matrix&, const psi::Psi<std::complex<double>>&, DensityMatrix<std::complex<double>, std::complex<double>>&);
}
