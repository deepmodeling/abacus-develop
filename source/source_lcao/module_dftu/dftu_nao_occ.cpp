#include "dftu_nao_occ.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_nao_folding.h"
#include "source_base/timer.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_reduce.h"
#include "source_estate/occ_matrix.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/module_symmetry/symm_rotation_k.h"

// cal_occ_mat_k / cal_occ_mat_gamma take Plus_U_Base& dftu directly and read all
// occupation-matrix state (occ/save arrays, lookup table, nspin/npol, and the
// occmat_ready flag) from dftu.occmat() and the Plus_U_Base accessors.

namespace
{
// (symmetry) lazily-built rotation machinery, shared across SCF iterations
// of one run: process-lifetime static since cal_occ_mat_k has no natural
// per-ion-step owning object to hang this off (unlike the dft_plus_u=1
// operator path, which owns its own copy).
ModuleSymmetry::Symmetry_rotation_k dftu_occ_symrot;
bool dftu_occ_symrot_built = false;

// Extract the real part (double -> identity, complex<double> -> .real()).
inline double real_part(const double& x) { return x; }
inline double real_part(const std::complex<double>& x) { return x.real(); }
// Extract the real part of the conjugate (double -> identity,
// complex<double> -> conj(x).real() == x.real()).
inline double real_part_conj(const double& x) { return x; }
inline double real_part_conj(const std::complex<double>& x)
{
    return std::conj(x).real();
}

/// @brief Accumulate one (iat, l, spin) channel of the occupation matrix
///        from the S*DM product srho. Reads npol and the corr_iwt lookup
///        directly from occmat so callers do not thread those scalars through.
template <typename T>
void acc_channel(OccupationMatrix& occmat,
                            const Parallel_Orbitals& pv,
                            const T* srho,
                            int iat,
                            int l,
                            int spin)
{
    const int npol = occmat.npol();
    ModuleBase::matrix& occ = occmat.mat(iat, l, spin);
    const int two_l_plus_one = 2 * l + 1;
    for (int m0 = 0; m0 < two_l_plus_one; m0++)
    {
        for (int ipol0 = 0; ipol0 < npol; ipol0++)
        {
            const int iwt0 = occmat.corr_iwt(iat, l, m0, ipol0);
            const int mu = pv.global2local_row(iwt0);
            const int mu_prime = pv.global2local_col(iwt0);

            for (int m1 = 0; m1 < two_l_plus_one; m1++)
            {
                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                {
                    const int iwt1 = occmat.corr_iwt(iat, l, m1, ipol1);
                    const int nu = pv.global2local_col(iwt1);
                    const int nu_prime = pv.global2local_row(iwt1);

                    const int irc = nu * pv.nrow + mu;
                    const int irc_prime = mu_prime * pv.nrow + nu_prime;

                    const int m0_all = m0 + ipol0 * two_l_plus_one;
                    const int m1_all = m1 + ipol1 * two_l_plus_one;

                    if ((nu >= 0) && (mu >= 0))
                    {
                        occ(m0_all, m1_all) += real_part(srho[irc]) / 4.0;
                    }

                    if ((nu_prime >= 0) && (mu_prime >= 0))
                    {
                        occ(m0_all, m1_all)
                            += real_part_conj(srho[irc_prime]) / 4.0;
                    }
                } // ipol1
            } // m1
        } // ipol0
    } // m0
}

/// @brief Walk the (it, ia, l) atom mesh for one k-point and accumulate
///        each qualifying channel of occmat from the S*DM product srho.
template <typename T>
void acc_for_ik(OccupationMatrix& occmat,
                           const UnitCell& ucell,
                           const Parallel_Orbitals& pv,
                           const T* srho,
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

                acc_channel(occmat, pv, srho, iat, l, spin);
            } // end l
        } // end ia
    } // end it
}

/// @brief accumulate one k-star member's rotated S*DM product into occmat,
///        redistributing the ibz k-point's full weight (already baked into
///        srho_ibz) across all kstar_size members via Symmetry_rotation's
///        built-in 1/kstar_size scaling (see restore_dm/rot_matrix_ao).
void acc_over_kstar(OccupationMatrix& occmat,
                               const UnitCell& ucell,
                               const Parallel_Orbitals& pv,
                               const K_Vectors& kv,
                               const std::vector<std::complex<double>>& srho_ibz,
                               const int ik_ibz,
                               const int spin,
                               const int nspin,
                               const std::vector<int>& l_channel)
{
    const int nsym = ucell.symm.nrotk;
    const size_t kstar_size = kv.kstars[ik_ibz].size();
    std::vector<std::complex<double>> sigma_y;
    for (const std::pair<const int, ModuleBase::Vector3<double>>& isym_kvd : kv.kstars[ik_ibz])
    {
        const int isym = isym_kvd.first;
        std::vector<std::complex<double>> srho_rot;
        if (isym < nsym)
        { // unitary space-group operation (isym==0 is the identity)
            srho_rot = dftu_occ_symrot.rot_matrix_ao(srho_ibz, ik_ibz, kstar_size, isym, pv);
        }
        else
        { // antiunitary element: TRS * (spatial operation), see restore_dm
            const int isym_M = ucell.symm.magnetic_nspin4 ? isym : (isym - nsym);
            if (nspin == 4)
            {
                if (sigma_y.empty()) { sigma_y = dftu_occ_symrot.set_sigma_y_2d(pv); }
                srho_rot = dftu_occ_symrot.trs_spin_rotate(
                    dftu_occ_symrot.rot_matrix_ao(srho_ibz, ik_ibz, kstar_size, isym_M, pv, false),
                    sigma_y, pv, 1.0);
            }
            else
            {
                srho_rot = dftu_occ_symrot.rot_matrix_ao(srho_ibz, ik_ibz, kstar_size, isym_M, pv, true);
            }
        }
        acc_for_ik(occmat, ucell, pv, srho_rot.data(), spin, l_channel);
    }
}

/// @brief MPI Allreduce each (iat, l) channel of occmat across all ranks
///        and symmetrize it (Hermitian average) per the nspin convention:
///        nspin=1 mirrors spin-0 into spin-1; nspin=2 symmetrizes each spin;
///        nspin=4 symmetrizes the single Pauli block. Reads nspin and npol
///        from occmat so callers do not thread them through.
void reduce_symm(OccupationMatrix& occmat,
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

                if (nspin == 1 || nspin == 4)
                {
                    ModuleBase::matrix& occ0 = occmat.mat(iat, l, 0);
                    // MPI Allreduce across ranks (in-place)
                    Parallel_Reduce::reduce_all(&occ0(0, 0),
                                                (2 * l + 1) * npol * (2 * l + 1) * npol);
                }
                else if (nspin == 2)
                {
                    ModuleBase::matrix& occ0 = occmat.mat(iat, l, 0);
                    // MPI Allreduce across ranks (in-place)
                    Parallel_Reduce::reduce_all(&occ0(0, 0),
                                                (2 * l + 1) * (2 * l + 1));

                    ModuleBase::matrix& occ1 = occmat.mat(iat, l, 1);
                    // MPI Allreduce across ranks (in-place)
                    Parallel_Reduce::reduce_all(&occ1(0, 0),
                                                (2 * l + 1) * (2 * l + 1));
                }

                switch (nspin)
                {
                case 1:
                {
                    ModuleBase::matrix& occ0 = occmat.mat(iat, l, 0);
                    occ0 += transpose(occ0);
                    occ0 *= 0.5;
                    occmat.mat(iat, l, 1) += occ0;
                    break;
                }

                case 2:
                    for (int is = 0; is < nspin; is++)
                    {
                        ModuleBase::matrix& occ_is = occmat.mat(iat, l, is);
                        occ_is += transpose(occ_is);
                    }
                    break;

                case 4:
                {
                    ModuleBase::matrix& occ0 = occmat.mat(iat, l, 0);
                    occ0 += transpose(occ0);
                    break;
                }

                default:
                    ModuleBase::WARNING_QUIT("DFTU_LCAO", "Not supported NSPIN parameter");
                }
            } // end l
        } // end ia
    } // end it
}

/// @brief Process one (it, ia, l, spin) block of the gamma-only
///        occupation matrix: accumulate from the real S*DM product srho,
///        MPI-Allreduce across ranks, then symmetrize per the nspin
///        convention. Reads nspin and npol from occmat so callers do not
///        thread them through.
void acc_channel_gamma(OccupationMatrix& occmat,
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

                // Calculate the local occupation number matrix
                acc_channel(occmat, pv, srho, iat, l, spin);
                ModuleBase::matrix& occ_is = occmat.mat(iat, l, spin);

                // MPI Allreduce across ranks (in-place)
                Parallel_Reduce::reduce_all(&occ_is(0, 0),
                                            (2 * l + 1) * npol * (2 * l + 1) * npol);

                // for the case spin independent calculation
                switch (nspin)
                {
                case 1:
                {
                    ModuleBase::matrix& occ0 = occmat.mat(iat, l, 0);
                    occ0 += transpose(occ0);
                    occ0 *= 0.5;
                    occmat.mat(iat, l, 1) += occ0;
                    break;
                }

                case 2:
                    occ_is += transpose(occ_is);
                    break;

                default:
                    ModuleBase::WARNING_QUIT("DFTU_LCAO", "Not supported NSPIN parameter");
                }
            } // L
        } // ia
    } // it
}

} // namespace

namespace DFTU_LCAO {

void cal_occ_mat_k(const Parallel_Orbitals* pv,
                   const UnitCell& ucell,
                   const std::vector<std::vector<std::complex<double>>>& dm_k,
                   const K_Vectors& kv,
                   const double& mixing_beta,
                   hamilt::Hamilt<std::complex<double>>* p_ham,
                   const bool gamma_only_local,
                   Plus_U_Base& dftu,
                   const std::string& ks_solver)
{
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

    // (symmetry) when crystal symmetry reduces the k-mesh, each ik below is only
    // the irreducible representative; build the AO rotation machinery once so
    // its k-star can be correctly re-expanded (see accumulate_occ_over_kstar).
    // Symmetry is analyzed once at the beginning and preserved by symmetrization.
    // Accordingly, symrot_, dftu_occ_symrot, and the cached Ms_ remain valid and 
    // do not need to be rebuilt each ionic step.
    const bool dftu_spacegroup_symmetry = (ModuleSymmetry::Symmetry::symm_flag == 1) && !kv.kstars.empty();
    if (dftu_spacegroup_symmetry && !dftu_occ_symrot_built)
    {
        const std::array<int, 3> period{ kv.nmp[0], kv.nmp[1], kv.nmp[2] };
        dftu_occ_symrot.find_irreducible_sector(ucell.symm, ucell.atoms, ucell.st,
            ModuleSymmetry::Symmetry_rotation_k::get_bvk_cells(period), period, ucell.lat);
        dftu_occ_symrot.cal_Ms(kv, ucell, *pv, nspin);
        dftu_occ_symrot_built = true;
    }

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
        if (dftu_spacegroup_symmetry)
        {
            // kv.kstars/Ms_ are sized per spin and indexed by GLOBAL ibz position
            // (kv.kstars.size() == nks_ibz); ik is local to this k-point pool, so
            // map it to the global k index first (kv.ik2iktot), then wrap into the
            // per-spin ibz range (mirrors RI_2D_Comm::split_m2D_ktoR_k's
            // "ik % ik_list.size()", but on the global index rather than the local one).
            const int ik_ibz = kv.ik2iktot[ik] % static_cast<int>(kv.kstars.size());
            acc_over_kstar(dftu.occmat(), ucell, *pv, kv, srho, ik_ibz, spin, nspin, l_channel);
        }
        else
        {
            acc_for_ik(dftu.occmat(), ucell, *pv, srho.data(), spin, l_channel);
        }
    } // ik

    // MPI Allreduce + symmetrize per (iat, l, n=0) channel across all ranks
    reduce_symm(dftu.occmat(), ucell, l_channel);

    if(dftu.has_occ_mixer() && dftu.is_occmat_ready())
    {
        dftu.occ_mixer().mix_plain(dftu.occmat(), mixing_beta);
    }

    dftu.set_occmat_ready();
    ModuleBase::timer::end("DFTU_LCAO", "cal_occ_mat_k");
    return;
}

void cal_occ_mat_gamma(const Parallel_Orbitals* pv,
                       const UnitCell &ucell,
                       const std::vector<std::vector<double>> &dm_gamma,
                       const double& mixing_beta,
                       hamilt::Hamilt<double>* p_ham,
                       Plus_U_Base& dftu)
{
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
    char transN = 'N';
    char transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0;
    const double beta = 0.0;

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
        acc_channel_gamma(dftu.occmat(), ucell, *pv, srho.data(), is, l_channel);
    } // is

    if(dftu.has_occ_mixer() && dftu.is_occmat_ready())
    {
        dftu.occ_mixer().mix_plain(dftu.occmat(), mixing_beta);
    }

    dftu.set_occmat_ready();
    ModuleBase::timer::end("DFTU_LCAO", "cal_occ_mat_gamma");
    return;
}

} // namespace DFTU_LCAO

//! dftu occupation matrix for multi-k case using dm(double)
template <>
void DFTU_LCAO::cal_occ_mat(const Parallel_Orbitals* pv,
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
    cal_occ_mat_k(pv, ucell, dm, kv, mixing_beta, p_ham, gamma_only_local, dftu, ks_solver);
}

//! dftu occupation matrix for gamma only using dm(double)
template <>
void DFTU_LCAO::cal_occ_mat(const Parallel_Orbitals* pv,
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
    cal_occ_mat_gamma(pv, ucell, dm, mixing_beta, p_ham, dftu);
}

template void DFTU_LCAO::cal_occ_mat(const Parallel_Orbitals*,
                                     const UnitCell&,
                                     const std::vector<std::vector<double>>&,
                                     const K_Vectors&,
                                     const double&,
                                     hamilt::Hamilt<double>*,
                                     Plus_U_Base&,
                                     const bool,
                                     const int,
                                     const std::string&);

template void DFTU_LCAO::cal_occ_mat(const Parallel_Orbitals*,
                                     const UnitCell&,
                                     const std::vector<std::vector<std::complex<double>>>&,
                                     const K_Vectors&,
                                     const double&,
                                     hamilt::Hamilt<std::complex<double>>*,
                                     Plus_U_Base&,
                                     const bool,
                                     const int,
                                     const std::string&);
