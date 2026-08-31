#include "dftu_lcao.h"
#include "dftu_hamilt.h"
#include "dftu_lcao_pots.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"


#ifdef __LCAO
namespace DFTU_LCAO {

void pot_uterm_complex(Plus_U& dftu,
                       const UnitCell& ucell,
                       const Parallel_Orbitals* pv,
                       const int ik,
                       std::complex<double>* pot_uterm,
                       const std::vector<int>& isk,
                       const std::complex<double>* sk,
                       const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_uterm_complex");
    if (!dftu.is_occ_mat_initialized())
    {
        return;
    }

    ModuleBase::timer::start("DFTU_LCAO", "pot_uterm_complex");

    int spin = isk[ik];

    const int nlocal = pv->get_global_row_size();
    ModuleBase::GlobalFunc::ZEROS(pot_uterm, pv->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> one(1.0, 0.0);
    const std::complex<double> half = 0.5;
    const std::complex<double> zero = 0.0;

    std::vector<std::complex<double>> pot_onsite(pv->nloc);
    DFTU_LCAO::pot_onsite_complex(dftu, ucell, pv, spin, true, &pot_onsite[0], npol);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), one_int, one_int, pv->desc,
            sk, one_int, one_int, pv->desc,
            zero,
            pot_uterm, one_int, one_int, pv->desc);
#endif

    for (int irc = 0; irc < pv->nloc; irc++)
    {
        pot_onsite[irc] = pot_uterm[irc];
    }

#ifdef __MPI
    ScalapackConnector::tranu(nlocal, nlocal,
            one,
            &pot_onsite[0], one_int, one_int, pv->desc,
            one,
            pot_uterm, one_int, one_int, pv->desc);
#endif

    ModuleBase::timer::end("DFTU_LCAO", "pot_uterm_complex");
    return;
}

void pot_uterm_real(Plus_U& dftu,
                    const UnitCell& ucell,
                    const Parallel_Orbitals* pv,
                    const int ik,
                    double* pot_uterm,
                    const std::vector<int>& isk,
                    const double* sk,
                    const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_uterm_real");
    if (!dftu.is_occ_mat_initialized())
    {
        return;
    }
    ModuleBase::timer::start("DFTU_LCAO", "pot_uterm_real");

    int spin = isk[ik];

    const int nlocal = pv->get_global_row_size();
    ModuleBase::GlobalFunc::ZEROS(pot_uterm, pv->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    int one_int = 1;
    double alpha = 1.0, beta = 0.0, half = 0.5, one = 1.0;

    std::vector<double> pot_onsite(pv->nloc);
    DFTU_LCAO::pot_onsite_real(dftu, ucell, pv, spin, true, &pot_onsite[0], npol);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), 1, 1, pv->desc,
            sk, 1, 1, pv->desc,
            beta,
            pot_uterm, 1, 1, pv->desc);
#endif

    for (int irc = 0; irc < pv->nloc; irc++)
        pot_onsite[irc] = pot_uterm[irc];

#ifdef __MPI
    pdtran_(&nlocal, &nlocal,
            &one,
            &pot_onsite[0], &one_int, &one_int, const_cast<int*>(pv->desc),
            &one,
            pot_uterm, &one_int, &one_int, const_cast<int*>(pv->desc));
#endif

    ModuleBase::timer::end("DFTU_LCAO", "pot_uterm_real");
    return;
}

} // namespace DFTU_LCAO

void Plus_U::cal_eff_pot_mat_R_double(const UnitCell& ucell, const Parallel_Orbitals* pv, const int ispin, double* SR, double* HR, const int npol)
{
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0, beta = 0.0, one = 1.0, half = 0.5;
    const int nlocal = pv->get_global_row_size();

    std::vector<double> pot_onsite(pv->nloc);
    DFTU_LCAO::pot_onsite_real(*this, ucell, pv, ispin, true, &pot_onsite[0], npol);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), 1, 1, pv->desc,
            SR, 1, 1, pv->desc,
            beta,
            HR, 1, 1, pv->desc);

    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            SR, 1, 1, pv->desc,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), 1, 1, pv->desc,
            one,
            HR, 1, 1, pv->desc);
#endif

    return;
}

void Plus_U::cal_eff_pot_mat_R_complex_double(const UnitCell& ucell, const Parallel_Orbitals* pv, const int ispin, std::complex<double>* SR, std::complex<double>* HR, const int npol)
{
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> zero = 0.0, one = 1.0, half = 0.5;
    const int nlocal = pv->get_global_row_size();

    std::vector<std::complex<double>> pot_onsite(pv->nloc);
    DFTU_LCAO::pot_onsite_complex(*this, ucell, pv, ispin, true, &pot_onsite[0], npol);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), one_int, one_int, pv->desc,
            SR, one_int, one_int, pv->desc,
            zero,
            HR, one_int, one_int, pv->desc);

    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            SR, one_int, one_int, pv->desc,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), one_int, one_int, pv->desc,
            one,
            HR, one_int, one_int, pv->desc);
#endif

    return;
}

#endif
