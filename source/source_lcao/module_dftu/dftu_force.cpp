#ifdef __LCAO
#include "dftu_force.h"
#include "dftu_folding.h"
#include "dftu_lcao.h"
#include "source_base/global_function.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

#include <complex>
#include <string>


namespace DFTU_LCAO {

void force_stress(Plus_U& dftu,
                  const UnitCell& ucell,
                  const Grid_Driver& gd,
                  std::vector<std::vector<double>>* dmk_d,
                  std::vector<std::vector<std::complex<double>>>* dmk_c,
                  const Parallel_Orbitals& pv,
                  ForceStressArrays& fsr,
                  ModuleBase::matrix& force_dftu,
                  ModuleBase::matrix& stress_dftu,
                  const K_Vectors& kv,
                  const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "force_stress");
    ModuleBase::timer::start("DFTU_LCAO", "force_stress");

    // Defensive null check: the legacy dft_plus_u==2 force/stress path
    // requires fsr.DSloc_x/y/z (gamma_only) or fsr.DSloc_Rx/Ry/Rz (multik)
    // and fsr.DH_r to be allocated and filled by the caller. If the caller
    // forgot to allocate them (as in force_stress_lcao.cpp where the local
    // fsr_dftu is created without allocation), we fail early with a clear
    // message instead of letting pdgemm_ dereference nullptr and crash.
    // See force_stress_lcao.cpp for the historical background.
    if (dftu.is_gamma_only_local())
    {
        if (dftu.is_cal_force()
            && (fsr.DSloc_x == nullptr || fsr.DSloc_y == nullptr || fsr.DSloc_z == nullptr))
        {
            ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress",
                "fsr.DSloc_x/y/z are nullptr in gamma_only path; the caller must allocate and fill them. "
                "See notes in source/source_lcao/force_stress_lcao.cpp.");
        }
        if (dftu.is_cal_stress()
            && (fsr.DSloc_x == nullptr || fsr.DSloc_y == nullptr || fsr.DSloc_z == nullptr
                || fsr.DH_r == nullptr))
        {
            ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress",
                "fsr.DSloc_x/y/z or fsr.DH_r is nullptr in gamma_only path; "
                "the caller must allocate and fill them. "
                "See notes in source/source_lcao/force_stress_lcao.cpp.");
        }
    }
    else
    {
        if (dftu.is_cal_force()
            && (fsr.DSloc_Rx == nullptr || fsr.DSloc_Ry == nullptr || fsr.DSloc_Rz == nullptr))
        {
            ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress",
                "fsr.DSloc_Rx/Ry/Rz are nullptr in multik path; the caller must allocate and fill them. "
                "See notes in source/source_lcao/force_stress_lcao.cpp.");
        }
        if (dftu.is_cal_stress()
            && (fsr.DSloc_Rx == nullptr || fsr.DSloc_Ry == nullptr || fsr.DSloc_Rz == nullptr
                || fsr.DH_r == nullptr))
        {
            ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress",
                "fsr.DSloc_Rx/Ry/Rz or fsr.DH_r is nullptr in multik path; "
                "the caller must allocate and fill them. "
                "See notes in source/source_lcao/force_stress_lcao.cpp.");
        }
    }

    // Layout invariant: the folded dSR/mat buffers are consumed by ScaLAPACK
    // GEMM through pv.desc (local column-major storage) and read back with
    // explicit ic * pv.nrow + ir indices. All ks_solvers accepted by INPUT
    // validation are column-major today; abort loudly instead of silently
    // producing wrong forces/stresses if that assumption ever changes.
    if ((dftu.is_cal_force() || dftu.is_cal_stress())
        && !ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(dftu.get_ks_solver()))
    {
        ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress",
            "non column-major ks_solver is not supported for DFT+U force/stress; "
            "the folded matrix layout assumption would be violated");
    }

    const int nlocal = dftu.get_nlocal();

    if (dftu.is_cal_force())
    {
        force_dftu.zero_out();
    }
    if (dftu.is_cal_stress())
    {
        stress_dftu.zero_out();
    }

    if (dftu.is_gamma_only_local())
    {
        const char transN = 'N';
        const char transT = 'T';
        const int one_int = 1;
        const double alpha = 1.0;
        const double beta = 0.0;

        std::vector<double> rho_pot_onsite(pv.nloc);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {

            const int spin = kv.isk[ik];

            double* pot_onsite = new double[pv.nloc];

            dftu.pot_onsite_real(spin, false, pot_onsite, npol);

#ifdef __MPI
            ScalapackConnector::gemm(transT, transN, nlocal, nlocal, nlocal,
                    alpha, (*dmk_d)[spin].data(), 1, 1,
                    pv.desc, pot_onsite, 1, 1,
                    pv.desc, beta, &rho_pot_onsite[0],
                    1, 1, pv.desc);
#endif

            delete[] pot_onsite;

            if (dftu.is_cal_force())
            {
                cal_force_gamma(dftu.get_nlocal(), dftu.get_npol(),
                                dftu.get_orbital_corr_vec(), dftu.get_iatlnmipol2iwt(),
                                ucell, &rho_pot_onsite[0], pv,
                                fsr.DSloc_x, fsr.DSloc_y, fsr.DSloc_z, force_dftu);
            }

            if (dftu.is_cal_stress())
            {
                cal_stress_gamma(dftu.get_nlocal(), dftu.get_npol(),
                                 dftu.get_ks_solver(), dftu.get_orb_cutoff(),
                                 ucell, pv, &gd,
                                 fsr.DSloc_x, fsr.DSloc_y, fsr.DSloc_z, fsr.DH_r,
                                 &rho_pot_onsite[0], stress_dftu);
            }
        } // ik
    }
    else
    {
        const char transN = 'N';
        const char transT = 'T';
        const int one_int = 1;
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);

        std::vector<std::complex<double>> rho_pot_onsite(pv.nloc);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            const int spin = kv.isk[ik];

            std::complex<double>* pot_onsite = new std::complex<double>[pv.nloc];

            dftu.pot_onsite_complex(spin, false, pot_onsite, npol);


#ifdef __MPI
            ScalapackConnector::gemm(transT, transN, nlocal, nlocal, nlocal,
                    alpha, (*dmk_c)[ik].data(), one_int, one_int,
                    pv.desc, pot_onsite, one_int, one_int, pv.desc, beta,
                    &rho_pot_onsite[0], one_int, one_int, pv.desc);
#endif

            delete[] pot_onsite;

            if (dftu.is_cal_force())
            {
                cal_force_k(dftu.get_nlocal(), dftu.get_npol(),
                            dftu.get_ks_solver(), dftu.get_orb_cutoff(),
                            dftu.get_orbital_corr_vec(), dftu.get_iatlnmipol2iwt(),
                            ucell, gd, fsr, pv, ik, &rho_pot_onsite[0], force_dftu, kv.kvec_d[ik]);
            }
            if (dftu.is_cal_stress())
            {
                cal_stress_k(dftu.get_nlocal(), dftu.get_npol(),
                             dftu.get_ks_solver(), dftu.get_orb_cutoff(),
                             ucell, gd, fsr, pv, ik, &rho_pot_onsite[0], stress_dftu, kv.kvec_d[ik]);
            }
        } // ik
    }

    if (dftu.is_cal_force())
    {
        Parallel_Reduce::reduce_pool(force_dftu.c, force_dftu.nr * force_dftu.nc);
    }

    if (dftu.is_cal_stress())
    {
        Parallel_Reduce::reduce_pool(stress_dftu.c, stress_dftu.nr * stress_dftu.nc);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i > j)
                    stress_dftu(i, j) = stress_dftu(j, i);
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                stress_dftu(i, j) *= ucell.lat0 / ucell.omega;
            }
        }
    }
    ModuleBase::timer::end("DFTU_LCAO", "force_stress");

    return;
}

void cal_force_k(const int nlocal,
                 const int npol,
                 const std::string& ks_solver,
                 const std::vector<double>& orb_cutoff,
                 const std::vector<int>& orbital_corr,
                 const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt,
                 const UnitCell& ucell,
                 const Grid_Driver& gd,
                 ForceStressArrays& fsr,
                 const Parallel_Orbitals& pv,
                 const int ik,
                 const std::complex<double>* rho_pot_onsite,
                 ModuleBase::matrix& force_dftu,
                 const ModuleBase::Vector3<double>& kvec_d)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_force_k");
    ModuleBase::timer::start("DFTU_LCAO", "cal_force_k");

    const char transN = 'N';
    const char transC = 'C';
    const int one_int = 1;
    const std::complex<double> zero(0.0, 0.0);
    const std::complex<double> one(1.0, 0.0);

    assert(nlocal>0);

    std::vector<std::complex<double>> dm_pot_onsite_dSm(pv.nloc);
    std::vector<std::complex<double>> dSm_k(pv.nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        DFTU_LCAO::folding_matrix_k(npol, ks_solver, orb_cutoff,
                                        ucell, gd, fsr, pv, ik, dim + 1, 0, &dSm_k[0], kvec_d);

#ifdef __MPI
        ScalapackConnector::gemm(transN,
                transC,
                nlocal,
                nlocal,
                nlocal,
                one,
                &dSm_k[0],
                one_int,
                one_int,
                pv.desc,
                rho_pot_onsite,
                one_int,
                one_int,
                pv.desc,
                zero,
                &dm_pot_onsite_dSm[0],
                one_int,
                one_int,
                pv.desc);
#endif

        for (int ir = 0; ir < pv.nrow; ir++)
        {
            const int iwt1 = pv.local2global_row(ir);
            const int iat1 = ucell.iwt2iat[iwt1];

            for (int ic = 0; ic < pv.ncol; ic++)
            {
                const int iwt2 = pv.local2global_col(ic);
                const int irc = ic * pv.nrow + ir;

                if (iwt1 == iwt2)
                    force_dftu(iat1, dim) += dm_pot_onsite_dSm[irc].real();

            } // end ic
        }     // end ir

#ifdef __MPI
        ScalapackConnector::gemm(transN,
                transN,
                nlocal,
                nlocal,
                nlocal,
                one,
                &dSm_k[0],
                one_int,
                one_int,
                pv.desc,
                rho_pot_onsite,
                one_int,
                one_int,
                pv.desc,
                zero,
                &dm_pot_onsite_dSm[0],
                one_int,
                one_int,
                pv.desc);
#endif

        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1)
                continue;
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it])
                        continue;
                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                            continue;

                        for (int m = 0; m < 2 * l + 1; m++)
                        {
                            for (int ipol = 0; ipol < npol; ipol++)
                            {
                                const int iwt = iatlnmipol2iwt[iat][l][n][m][ipol];
                                const int mu = pv.global2local_row(iwt);
                                const int nu = pv.global2local_col(iwt);
                                if (mu < 0 || nu < 0)
                                    continue;

                                force_dftu(iat, dim) += dm_pot_onsite_dSm[nu * pv.nrow + mu].real();
                            }
                        } //
                    }     // n
                }         // l
            }             // ia
        }                 // it
    }                     // end dim
    ModuleBase::timer::end("DFTU_LCAO", "cal_force_k");

    return;
}

void cal_stress_k(const int nlocal,
                  const int npol,
                  const std::string& ks_solver,
                  const std::vector<double>& orb_cutoff,
                  const UnitCell& ucell,
                  const Grid_Driver& gd,
                  ForceStressArrays& fsr,
                  const Parallel_Orbitals& pv,
                  const int ik,
                  const std::complex<double>* rho_pot_onsite,
                  ModuleBase::matrix& stress_dftu,
                  const ModuleBase::Vector3<double>& kvec_d)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_stress_k");
    ModuleBase::timer::start("DFTU_LCAO", "cal_stress_k");

    const char transN = 'N';
    const int one_int = 1;
    const std::complex<double> minus_half(-0.5, 0.0);
    const std::complex<double> zero(0.0, 0.0);
    const std::complex<double> one(1.0, 0.0);

    std::vector<std::complex<double>> dm_pot_onsite_sover(pv.nloc);
    std::vector<std::complex<double>> dSR_k(pv.nloc);

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            DFTU_LCAO::folding_matrix_k(npol, ks_solver, orb_cutoff,
                                            ucell, gd, fsr, pv, ik, dim1 + 4, dim2, &dSR_k[0], kvec_d);

#ifdef __MPI
            ScalapackConnector::gemm(transN,
                    transN,
                    nlocal,
                    nlocal,
                    nlocal,
                    minus_half,
                    rho_pot_onsite,
                    one_int,
                    one_int,
                    pv.desc,
                    &dSR_k[0],
                    one_int,
                    one_int,
                    pv.desc,
                    zero,
                    &dm_pot_onsite_sover[0],
                    one_int,
                    one_int,
                    pv.desc);
#endif

            for (int ir = 0; ir < pv.nrow; ir++)
            {
                const int iwt1 = pv.local2global_row(ir);
                for (int ic = 0; ic < pv.ncol; ic++)
                {
                    const int iwt2 = pv.local2global_col(ic);
                    const int irc = ic * pv.nrow + ir;

                    if (iwt1 == iwt2)
                        stress_dftu(dim1, dim2) += 2.0 * dm_pot_onsite_sover[irc].real();
                } // end ic
            }     // end ir

        } // end dim2
    }     // end dim1
    ModuleBase::timer::end("DFTU_LCAO", "cal_stress_k");

    return;
}

void cal_force_gamma(const int nlocal,
                     const int npol,
                     const std::vector<int>& orbital_corr,
                     const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt,
                     const UnitCell& ucell,
                     const double* rho_pot_onsite,
                     const Parallel_Orbitals& pv,
                     double* dsloc_x,
                     double* dsloc_y,
                     double* dsloc_z,
                     ModuleBase::matrix& force_dftu)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_force_gamma");
    ModuleBase::timer::start("DFTU_LCAO", "cal_force_gamma");
    const char transN = 'N';
    const char transT = 'T';
    const int one_int = 1;
    const double one = 1.0;
    const double zero = 0.0;
    const double minus_one = -1.0;
    assert(nlocal>0);

    std::vector<double> dm_pot_onsite_dSm(pv.nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        double* tmp_ptr = nullptr;
        if (dim == 0)
        {
            tmp_ptr = dsloc_x;
        }
        else if (dim == 1)
        {
            tmp_ptr = dsloc_y;
        }
        else if (dim == 2)
        {
            tmp_ptr = dsloc_z;
        }

#ifdef __MPI
        ScalapackConnector::gemm(transN,
                transT,
                nlocal,
                nlocal,
                nlocal,
                one,
                tmp_ptr,
                1,
                1,
                pv.desc,
                rho_pot_onsite,
                1,
                1,
                pv.desc,
                zero,
                &dm_pot_onsite_dSm[0],
                1,
                1,
                pv.desc);
#endif

        for (int ir = 0; ir < pv.nrow; ir++)
        {
            const int iwt1 = pv.local2global_row(ir);
            const int iat1 = ucell.iwt2iat[iwt1];

            for (int ic = 0; ic < pv.ncol; ic++)
            {
                const int iwt2 = pv.local2global_col(ic);
                const int irc = ic * pv.nrow + ir;

                if (iwt1 == iwt2)
                    force_dftu(iat1, dim) += dm_pot_onsite_dSm[irc];

            } // end ic
        }     // end ir

#ifdef __MPI
        ScalapackConnector::gemm(transN,
                transT,
                nlocal,
                nlocal,
                nlocal,
                one,
                tmp_ptr,
                1,
                1,
                pv.desc,
                rho_pot_onsite,
                1,
                1,
                pv.desc,
                zero,
                &dm_pot_onsite_dSm[0],
                1,
                1,
                pv.desc);
#endif

        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1)
                continue;
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it])
                        continue;

                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                            continue;

                        // Calculate the local occupation number matrix
                        for (int m = 0; m < 2 * l + 1; m++)
                        {
                            for (int ipol = 0; ipol < npol; ipol++)
                            {
                                const int iwt = iatlnmipol2iwt[iat][l][n][m][ipol];
                                const int mu = pv.global2local_row(iwt);
                                const int nu = pv.global2local_col(iwt);
                                if (mu < 0 || nu < 0)
                                    continue;

                                force_dftu(iat, dim) += dm_pot_onsite_dSm[nu * pv.nrow + mu];
                            }
                        } //
                    }     // n
                }         // l
            }             // ia
        }                 // it

    } // end dim
    ModuleBase::timer::end("DFTU_LCAO", "cal_force_gamma");

    return;
}

void cal_stress_gamma(const int nlocal,
                      const int npol,
                      const std::string& ks_solver,
                      const std::vector<double>& orb_cutoff,
                      const UnitCell& ucell,
                      const Parallel_Orbitals& pv,
                      const Grid_Driver* gd,
                      double* dsloc_x,
                      double* dsloc_y,
                      double* dsloc_z,
                      double* dh_r,
                      const double* rho_pot_onsite,
                      ModuleBase::matrix& stress_dftu)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_stress_gamma");
    ModuleBase::timer::start("DFTU_LCAO", "cal_stress_gamma");

    const char transN = 'N';
    const int one_int = 1;
    const double zero = 0.0;
    const double minus_half = -0.5;
    const double one = 1.0;

    std::vector<double> dSR_gamma(pv.nloc);
    std::vector<double> dm_pot_onsite_sover(pv.nloc);

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            DFTU_LCAO::fold_dSR_gamma(npol, ks_solver, orb_cutoff,
                                         ucell, pv, gd, dsloc_x, dsloc_y, dsloc_z, dh_r, dim1, dim2, &dSR_gamma[0]);

#ifdef __MPI
            ScalapackConnector::gemm(transN,
                    transN,
                    nlocal,
                    nlocal,
                    nlocal,
                    minus_half,
                    rho_pot_onsite,
                    1,
                    1,
                    pv.desc,
                    &dSR_gamma[0],
                    1,
                    1,
                    pv.desc,
                    zero,
                    &dm_pot_onsite_sover[0],
                    1,
                    1,
                    pv.desc);
#endif

            for (int ir = 0; ir < pv.nrow; ir++)
            {
                const int iwt1 = pv.local2global_row(ir);

                for (int ic = 0; ic < pv.ncol; ic++)
                {
                    const int iwt2 = pv.local2global_col(ic);
                    const int irc = ic * pv.nrow + ir;

                    if (iwt1 == iwt2)
                        stress_dftu(dim1, dim2) += 2.0 * dm_pot_onsite_sover[irc];
                } // end ic
            }     // end ir

        } // end dim2
    }     // end dim1
    ModuleBase::timer::end("DFTU_LCAO", "cal_stress_gamma");
    return;
}

} // namespace DFTU_LCAO

#endif
