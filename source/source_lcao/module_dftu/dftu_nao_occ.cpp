#include "dftu_nao_occ.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_nao_folding.h"
#include "source_base/timer.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_reduce.h"
#include "source_estate/occ_matrix.h"
#include "source_lcao/hamilt_lcao.h"

// cal_occ_mat_k / cal_occ_mat_gamma take Plus_U_Base& dftu directly and read all
// occupation-matrix state (occ/save arrays, lookup table, nspin/npol, and the
// occmat_ready flag) from dftu.occmat() and the Plus_U_Base accessors.


void DFTU_LCAO::cal_occ_mat_k(const Parallel_Orbitals* pv,
                         const UnitCell& ucell,
                         const std::vector<std::vector<std::complex<double>>>& dm_k,
                         const K_Vectors& kv,
                         const double& mixing_beta,
                         hamilt::Hamilt<std::complex<double>>* p_ham,
                         const bool gamma_only_local,
                         Plus_U_Base& dftu,
                         const std::string& ks_solver)
{
    if (pv == nullptr || p_ham == nullptr)
    {
        return;
    }
    ModuleBase::TITLE("DFTU_LCAO", "cal_occ_mat_k");
    ModuleBase::timer::start("DFTU_LCAO", "cal_occ_mat_k");

    const int nspin = dftu.occmat().nspin();
    const int nlocal = pv->get_global_row_size();
    const std::vector<int>& l_channel = dftu.get_l_channel_vec();

    // copy occ_mat to occ_mat_save, then zero occ_mat
    dftu.occmat().copy_to_save(ucell, l_channel);
    dftu.occmat().zero(ucell, l_channel);

    //=================Part 1======================
    // call SCALAPACK routine to calculate the product of the S and density matrix
    const char transN = 'N';
    const char transT = 'T';
    const int one_int = 1;
    const std::complex<double> beta(0.0,0.0), alpha(1.0,0.0);

    std::vector<std::complex<double>> srho(pv->nloc);

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        // srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
        DFTU_LCAO::folding_matrix_k_new(ks_solver, gamma_only_local, nspin, ik, p_ham);

        std::complex<double>* s_k_pointer = nullptr;

        if(nspin != 4)
        {
            s_k_pointer = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)->getSk();
        }
        else
        {
            s_k_pointer = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->getSk();
        }

#ifdef __MPI
        ScalapackConnector::gemm(transN,
            transT,
            nlocal,
            nlocal,
            nlocal,
            alpha,
            s_k_pointer,
            one_int,
            one_int,
            &pv->desc[0],
            dm_k[ik].data(),
            one_int,
            one_int,
            &pv->desc[0],
            beta,
            srho.data(),
            one_int,
            one_int,
            &pv->desc[0]);
#endif

        const int spin = kv.isk[ik];
        // Walk (it, ia, l, n=0) and accumulate each qualifying channel
        accumulate_occ_k_for_ik(dftu.occmat(), ucell, *pv, srho.data(), spin, l_channel);
    } // ik

    // MPI Allreduce + symmetrize per (iat, l, n=0) channel across all ranks
    reduce_and_symmetrize_occ_k(dftu.occmat(), ucell, l_channel);

    if(dftu.has_occ_mixer() && dftu.is_occmat_ready())
    {
        dftu.occ_mixer().mix_plain(dftu.occmat(), mixing_beta);
    }

    dftu.set_occmat_ready();
    ModuleBase::timer::end("DFTU_LCAO", "cal_occ_mat_k");
    return;
}

void DFTU_LCAO::cal_occ_mat_gamma(const Parallel_Orbitals* pv,
                             const UnitCell &ucell,
                             const std::vector<std::vector<double>> &dm_gamma,
                             const double& mixing_beta,
                             hamilt::Hamilt<double>* p_ham,
                             Plus_U_Base& dftu)
{
    if (pv == nullptr || p_ham == nullptr)
    {
        return;
    }
    ModuleBase::TITLE("DFTU_LCAO", "cal_occ_mat_gamma");
    ModuleBase::timer::start("DFTU_LCAO", "cal_occ_mat_gamma");

    const int nspin = dftu.occmat().nspin();
    const int nlocal = pv->get_global_row_size();
    const std::vector<int>& l_channel = dftu.get_l_channel_vec();

    // copy occ_mat to occ_mat_save, then zero occ_mat
    dftu.occmat().copy_to_save(ucell, l_channel);
    dftu.occmat().zero(ucell, l_channel);

    //=================Part 1======================
    // call PBLAS routine to calculate the product of the S and density matrix
    char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0, beta = 0.0;

    std::vector<double> srho(pv->nloc);
    for (int is = 0; is < nspin; is++)
    {
        double* s_gamma_pointer = dynamic_cast<hamilt::HamiltLCAO<double, double>*>(p_ham)->getSk();

#ifdef __MPI
        ScalapackConnector::gemm(transN,
            transT,
            nlocal,
            nlocal,
            nlocal,
            alpha,
            s_gamma_pointer,
            one_int,
            one_int,
            &pv->desc[0],
            dm_gamma[is].data(),
            //dm_gamma[is].c,
            one_int,
            one_int,
            &pv->desc[0],
            beta,
            srho.data(),
            one_int,
            one_int,
            &pv->desc[0]);
#endif

        // Per (it, ia, l, n=0, spin) block: accumulate + Allreduce + symmetrize
        process_occ_channel_gamma(dftu.occmat(), ucell, *pv, srho.data(), is, l_channel);
    } // is

    if(dftu.has_occ_mixer() && dftu.is_occmat_ready())
    {
        dftu.occ_mixer().mix_plain(dftu.occmat(), mixing_beta);
    }

    dftu.set_occmat_ready();
    ModuleBase::timer::end("DFTU_LCAO", "cal_occ_mat_gamma");
    return;
}

namespace DFTU_LCAO {

/// @brief Accumulate one (iat, l, n, spin) channel of the occupation matrix
///        from the complex S*DM product srho for the multi-k case. Reads npol
///        and the iatlnmipol2iwt lookup directly from occmat so callers do
///        not need to thread those scalars through.
void accumulate_occ_channel_k(OccupationMatrix& occmat,
                              const Parallel_Orbitals& pv,
                              const std::complex<double>* srho,
                              int iat,
                              int l,
                              int n,
                              int spin)
{
    const int npol = occmat.npol();
    const auto& iatlnmipol2iwt = occmat.iatlnmipol2iwt();
    ModuleBase::matrix& occ = occmat.mat(iat, l, n, spin);
    const int two_l_plus_one = 2 * l + 1;
    for (int m0 = 0; m0 < two_l_plus_one; m0++)
    {
        for (int ipol0 = 0; ipol0 < npol; ipol0++)
        {
            const int iwt0 = iatlnmipol2iwt[iat][l][n][m0][ipol0];
            const int mu = pv.global2local_row(iwt0);
            const int mu_prime = pv.global2local_col(iwt0);

            for (int m1 = 0; m1 < two_l_plus_one; m1++)
            {
                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                {
                    const int iwt1 = iatlnmipol2iwt[iat][l][n][m1][ipol1];
                    const int nu = pv.global2local_col(iwt1);
                    const int nu_prime = pv.global2local_row(iwt1);

                    const int irc = nu * pv.nrow + mu;
                    const int irc_prime = mu_prime * pv.nrow + nu_prime;

                    const int m0_all = m0 + ipol0 * two_l_plus_one;
                    const int m1_all = m1 + ipol1 * two_l_plus_one;

                    if ((nu >= 0) && (mu >= 0))
                    {
                        occ(m0_all, m1_all) += (srho[irc]).real() / 4.0;
                    }

                    if ((nu_prime >= 0) && (mu_prime >= 0))
                    {
                        occ(m0_all, m1_all)
                            += (std::conj(srho[irc_prime])).real() / 4.0;
                    }
                } // ipol1
            } // m1
        } // ipol0
    } // m0
}

/// @brief Accumulate one (iat, l, n, spin) channel of the occupation matrix
///        from the real S*DM product srho for the gamma-only case. Reads npol
///        and the iatlnmipol2iwt lookup directly from occmat so callers do
///        not need to thread those scalars through. Uses the combined
///        (m0_all, m1_all) channel index consistently with the multi-k path.
void accumulate_occ_channel_gamma(OccupationMatrix& occmat,
                                  const Parallel_Orbitals& pv,
                                  const double* srho,
                                  int iat,
                                  int l,
                                  int n,
                                  int spin)
{
    const int npol = occmat.npol();
    const auto& iatlnmipol2iwt = occmat.iatlnmipol2iwt();
    ModuleBase::matrix& occ_is = occmat.mat(iat, l, n, spin);
    const int two_l_plus_one = 2 * l + 1;
    for (int m0 = 0; m0 < two_l_plus_one; m0++)
    {
        for (int ipol0 = 0; ipol0 < npol; ipol0++)
        {
            const int iwt0 = iatlnmipol2iwt[iat][l][n][m0][ipol0];
            const int mu = pv.global2local_row(iwt0);
            const int mu_prime = pv.global2local_col(iwt0);

            for (int m1 = 0; m1 < two_l_plus_one; m1++)
            {
                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                {
                    const int iwt1 = iatlnmipol2iwt[iat][l][n][m1][ipol1];
                    const int nu = pv.global2local_col(iwt1);
                    const int nu_prime = pv.global2local_row(iwt1);

                    const int irc = nu * pv.nrow + mu;
                    const int irc_prime = mu_prime * pv.nrow + nu_prime;

                    const int m0_all = m0 + ipol0 * two_l_plus_one;
                    const int m1_all = m1 + ipol1 * two_l_plus_one;

                    if ((nu >= 0) && (mu >= 0))
                    {
                        occ_is(m0_all, m1_all) += srho[irc] / 4.0;
                    }

                    if ((nu_prime >= 0) && (mu_prime >= 0))
                    {
                        occ_is(m0_all, m1_all) += srho[irc_prime] / 4.0;
                    }
                } // ipol1
            } // m1
        } // ipol0
    } // m0
}

/// @brief MPI Allreduce each (iat, l, n=0) channel of occmat across all ranks
///        and symmetrize it (Hermitian average) per the nspin convention:
///        nspin=1 mirrors spin-0 into spin-1; nspin=2 symmetrizes each spin;
///        nspin=4 symmetrizes the single Pauli block. Reads nspin and npol
///        from occmat so callers do not thread them through.
void reduce_and_symmetrize_occ_k(OccupationMatrix& occmat,
                                 const UnitCell& ucell,
                                 const std::vector<int>& l_channel)
{
    const int nspin = occmat.nspin();
    const int npol = occmat.npol();
    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;
        const int LC = l_channel[it];

        if (LC == -1)
        {
            continue;
        }

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);

            for (int l = 0; l < NL; l++)
            {
                if (l != l_channel[it])
                {
                    continue;
                }

                const int N = ucell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    // if(!Yukawa && n!=0) continue;
                    if (n != 0)
                    {
                        continue;
                    }
                    // set the local occupation mumber matrix of spin up and down zeros

                    if (nspin == 1 || nspin == 4)
                    {
                        ModuleBase::matrix& occ0 = occmat.mat(iat, l, n, 0);
                        // MPI Allreduce across ranks (in-place)
                        Parallel_Reduce::reduce_all(&occ0(0, 0),
                                                    (2 * l + 1) * npol * (2 * l + 1) * npol);
                    }
                    else if (nspin == 2)
                    {
                        ModuleBase::matrix& occ0 = occmat.mat(iat, l, n, 0);
                        // MPI Allreduce across ranks (in-place)
                        Parallel_Reduce::reduce_all(&occ0(0, 0),
                                                    (2 * l + 1) * (2 * l + 1));

                        ModuleBase::matrix& occ1 = occmat.mat(iat, l, n, 1);
                        // MPI Allreduce across ranks (in-place)
                        Parallel_Reduce::reduce_all(&occ1(0, 0),
                                                    (2 * l + 1) * (2 * l + 1));
                    }

                    switch (nspin)
                    {
                    case 1:
                    {
                        ModuleBase::matrix& occ0 = occmat.mat(iat, l, n, 0);
                        occ0 += transpose(occ0);
                        occ0 *= 0.5;
                        occmat.mat(iat, l, n, 1) += occ0;
                        break;
                    }

                    case 2:
                        for (int is = 0; is < nspin; is++)
                        {
                            ModuleBase::matrix& occ_is = occmat.mat(iat, l, n, is);
                            occ_is += transpose(occ_is);
                        }
                        break;

                    case 4:
                    {
                        ModuleBase::matrix& occ0 = occmat.mat(iat, l, n, 0);
                        occ0 += transpose(occ0);
                        break;
                    }

                    default:
                        std::cout << "Not supported NSPIN parameter" << std::endl;
                        exit(0);
                    }
                } // end n
            } // end l
        } // end ia
    } // end it
}

/// @brief Walk the (it, ia, l, n=0) atom mesh for one k-point and accumulate
///        each qualifying channel of occmat from the complex S*DM product
///        srho. Reads npol and the iatlnmipol2iwt lookup from occmat so
///        callers do not thread them through.
void accumulate_occ_k_for_ik(OccupationMatrix& occmat,
                             const UnitCell& ucell,
                             const Parallel_Orbitals& pv,
                             const std::complex<double>* srho,
                             int spin,
                             const std::vector<int>& l_channel)
{
    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;
        const int LC = l_channel[it];

        if (LC == -1)
        {
            continue;
        }

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);

            for (int l = 0; l < NL; l++)
            {
                if (l != l_channel[it])
                {
                    continue;
                }

                const int N = ucell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    // if(!Yukawa && n!=0) continue;
                    if (n != 0)
                    {
                        continue;
                    }

                    // Calculate the local occupation number matrix
                    accumulate_occ_channel_k(occmat, pv, srho, iat, l, n, spin);
                } // end n
            } // end l
        } // end ia
    } // end it
}

/// @brief Process one (it, ia, l, n=0, spin) block of the gamma-only
///        occupation matrix: accumulate from the real S*DM product srho,
///        MPI-Allreduce across ranks, then symmetrize per the nspin
///        convention. Reads nspin and npol from occmat so callers do not
///        thread them through.
void process_occ_channel_gamma(OccupationMatrix& occmat,
                               const UnitCell& ucell,
                               const Parallel_Orbitals& pv,
                               const double* srho,
                               int spin,
                               const std::vector<int>& l_channel)
{
    const int nspin = occmat.nspin();
    const int npol = occmat.npol();
    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;
        const int LC = l_channel[it];

        if (LC == -1)
        {
            continue;
        }
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);

            for (int l = 0; l < NL; l++)
            {
                if (l != l_channel[it])
                {
                    continue;
                }

                const int N = ucell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    // Calculate the local occupation number matrix
                    accumulate_occ_channel_gamma(occmat, pv, srho, iat, l, n, spin);
                    ModuleBase::matrix& occ_is = occmat.mat(iat, l, n, spin);

                    // MPI Allreduce across ranks (in-place)
                    Parallel_Reduce::reduce_all(&occ_is(0, 0),
                                                (2 * l + 1) * npol * (2 * l + 1) * npol);

                    // for the case spin independent calculation
                    switch (nspin)
                    {
                    case 1:
                    {
                        ModuleBase::matrix& occ0 = occmat.mat(iat, l, n, 0);
                        occ0 += transpose(occ0);
                        occ0 *= 0.5;
                        occmat.mat(iat, l, n, 1) += occ0;
                        break;
                    }

                    case 2:
                        occ_is += transpose(occ_is);
                        break;

                    default:
                        std::cout << "Not supported NSPIN parameter" << std::endl;
                        exit(0);
                    }

                } // end for(n)
            } // L
        } // ia
    } // it
}

//! dftu occupation matrix for gamma only using dm(double)
template <>
void cal_occ_mat(const Parallel_Orbitals* pv,
                 const UnitCell& ucell,
                 const std::vector<std::vector<double>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<double>* p_ham,
                 Plus_U_Base& dftu,
                 const bool gamma_only_local,
                 const int nspin,
                 const std::string& ks_solver)
{
    DFTU_LCAO::cal_occ_mat_gamma(pv, ucell, dm, mixing_beta, p_ham, dftu);
}

//! dftu occupation matrix for multiple k-points using dm(complex)
template <>
void cal_occ_mat(const Parallel_Orbitals* pv,
                 const UnitCell& ucell,
                 const std::vector<std::vector<std::complex<double>>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<std::complex<double>>* p_ham,
                 Plus_U_Base& dftu,
                 const bool gamma_only_local,
                 const int nspin,
                 const std::string& ks_solver)
{
    DFTU_LCAO::cal_occ_mat_k(pv, ucell, dm, kv, mixing_beta, p_ham, gamma_only_local, dftu, ks_solver);
}

} // namespace DFTU_LCAO
