#include "dftu_nao_fs_k.h"
#include "dftu_nao_folding.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_nao_pots.h"
#include "source_lcao/force_stress_arrays.h"
#include "source_base/global_function.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"

#include <complex>
#include <string>


namespace DFTU_LCAO {

namespace
{

/// @brief Add the real part of diagonal local-block entries to one force component.
///
/// Sums dm(ir, ic) over local block pairs whose global orbital indices
/// coincide, attributing each entry to the atom owning the orbital along
/// Cartesian component dim.
template <typename T>
void accumulate_diag_force(const Parallel_Orbitals& pv,
                           const UnitCell& ucell,
                           const T* dm,
                           const int dim,
                           ModuleBase::matrix& force_dftu)
{
    assert(dm != nullptr);
    assert(dim >= 0 && dim < 3);
    for (int ir = 0; ir < pv.nrow; ir++)
    {
        const int iwt1 = pv.local2global_row(ir);
        const int iat1 = ucell.iwt2iat[iwt1];
        for (int ic = 0; ic < pv.ncol; ic++)
        {
            if (pv.local2global_col(ic) == iwt1)
            {
                force_dftu(iat1, dim) += std::real(dm[ic * pv.nrow + ir]);
            }
        }
    }
}

/// @brief Add the real part of diagonal local-block entries to one stress pair.
template <typename T>
void accumulate_diag_stress(const Parallel_Orbitals& pv,
                            const T* dm,
                            const int dim1,
                            const int dim2,
                            const double factor,
                            ModuleBase::matrix& stress_dftu)
{
    assert(dm != nullptr);
    assert(dim1 >= 0 && dim1 < 3);
    assert(dim2 >= 0 && dim2 < 3);
    for (int ir = 0; ir < pv.nrow; ir++)
    {
        const int iwt1 = pv.local2global_row(ir);
        for (int ic = 0; ic < pv.ncol; ic++)
        {
            if (pv.local2global_col(ic) == iwt1)
            {
                stress_dftu(dim1, dim2) += factor * std::real(dm[ic * pv.nrow + ir]);
            }
        }
    }
}

/// @brief Whether atom type it carries a usable correlated channel.
///
/// Mirrors the silent skips of the original nest: no channel configured
/// (l_channel == -1), channel l beyond the atom's angular-momentum range,
/// or no n = 0 projector for that channel.
bool has_valid_correlated_channel(const UnitCell& ucell,
                                  const std::vector<int>& l_channel,
                                  const int it)
{
    const int lc = l_channel[it];
    return lc != -1 && lc < ucell.atoms[it].nwl + 1
           && ucell.atoms[it].l_nchi[lc] >= 1;
}

/// @brief Add the onsite (correlated-orbital) diagonal contribution to one force component.
///
/// For each type with a correlated channel, visits every atom and every
/// (m, spinor) orbital of the channel at the n = 0 projector. Equivalent
/// to the original it/ia/l/n/m/ipol nest, in which the l loop only ever
/// ran for l == l_channel and the n loop only for n == 0; the explicit
/// range guards preserve the silent skip when the channel is out of
/// range or the atom type has no n = 0 projector.
template <typename T>
void accumulate_onsite_force(Plus_U_Base& dftu,
                             const Parallel_Orbitals& pv,
                             const UnitCell& ucell,
                             const int npol,
                             const T* dm,
                             const int dim,
                             ModuleBase::matrix& force_dftu)
{
    assert(dm != nullptr);
    assert(dim >= 0 && dim < 3);
    assert(npol == 1 || npol == 2);
    const std::vector<int>& l_channel = dftu.get_l_channel_vec();
    const auto& iatlnmipol2iwt = dftu.occmat().iatlnmipol2iwt();
    for (int it = 0; it < ucell.ntype; it++)
    {
        if (!has_valid_correlated_channel(ucell, l_channel, it))
        {
            continue;
        }
        const int lc = l_channel[it];
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);
            for (int m = 0; m < 2 * lc + 1; m++)
            {
                for (int ipol = 0; ipol < npol; ipol++)
                {
                    const int iwt = iatlnmipol2iwt[iat][lc][0][m][ipol];
                    const int mu = pv.global2local_row(iwt);
                    const int nu = pv.global2local_col(iwt);
                    if (mu < 0 || nu < 0)
                    {
                        continue;
                    }
                    force_dftu(iat, dim) += std::real(dm[nu * pv.nrow + mu]);
                }
            }
        }
    }
}

} // namespace

void cal_force_k(const DftuFsEnv& env,
                 const int ik,
                 const ModuleBase::Vector3<double>& kvec_d,
                 const std::complex<double>* rho_pot_onsite,
                 ModuleBase::matrix& force_dftu)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_force_k");
    ModuleBase::timer::start("DFTU_LCAO", "cal_force_k");

    const Parallel_Orbitals& pv = env.pv();
    const UnitCell& ucell = env.ucell();
    const Grid_Driver& gd = env.gd();
    ForceStressArrays& fsr = env.fsr();
    const int npol = env.npol();
    const std::string& ks_solver = env.ks_solver();
    const std::vector<double>& orb_cutoff = env.orb_cutoff();
    const int nlocal = pv.get_global_row_size();

    // shared folding context: read-only params bundled for folding_matrix_k
    DFTU_LCAO::FoldingCtx fold_ctx{npol, ks_solver, orb_cutoff, &ucell, &pv, &gd};

    const char transN = 'N';
    const char transC = 'C';
    const int one_int = 1;
    const std::complex<double> zero(0.0, 0.0);
    const std::complex<double> one(1.0, 0.0);

    assert(nlocal > 0);

    std::vector<std::complex<double>> dm_pot_onsite_dSm(pv.nloc);
    std::vector<std::complex<double>> dSm_k(pv.nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        DFTU_LCAO::folding_matrix_k(fold_ctx, fsr, ik, dim + 1, 0, &dSm_k[0], kvec_d);

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

        accumulate_diag_force(pv, ucell, dm_pot_onsite_dSm.data(), dim, force_dftu);

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

        accumulate_onsite_force(env.dftu(), pv, ucell, npol,
                                dm_pot_onsite_dSm.data(), dim, force_dftu);
    }                     // end dim
    ModuleBase::timer::end("DFTU_LCAO", "cal_force_k");
}

void cal_stress_k(const DftuFsEnv& env,
                  const int ik,
                  const ModuleBase::Vector3<double>& kvec_d,
                  const std::complex<double>* rho_pot_onsite,
                  ModuleBase::matrix& stress_dftu)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_stress_k");
    ModuleBase::timer::start("DFTU_LCAO", "cal_stress_k");

    const Parallel_Orbitals& pv = env.pv();
    const UnitCell& ucell = env.ucell();
    const Grid_Driver& gd = env.gd();
    ForceStressArrays& fsr = env.fsr();
    const int npol = env.npol();
    const std::string& ks_solver = env.ks_solver();
    const std::vector<double>& orb_cutoff = env.orb_cutoff();
    const int nlocal = pv.get_global_row_size();

    // shared folding context: read-only params bundled for folding_matrix_k
    DFTU_LCAO::FoldingCtx fold_ctx{npol, ks_solver, orb_cutoff, &ucell, &pv, &gd};

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
            DFTU_LCAO::folding_matrix_k(fold_ctx, fsr, ik, dim1 + 4, dim2, &dSR_k[0], kvec_d);

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

            accumulate_diag_stress(pv, dm_pot_onsite_sover.data(), dim1, dim2, 2.0, stress_dftu);

        } // end dim2
    }     // end dim1
    ModuleBase::timer::end("DFTU_LCAO", "cal_stress_k");
}

void cal_force_gamma(const DftuFsEnv& env,
                     const double* rho_pot_onsite,
                     ModuleBase::matrix& force_dftu)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_force_gamma");
    ModuleBase::timer::start("DFTU_LCAO", "cal_force_gamma");

    const Parallel_Orbitals& pv = env.pv();
    const UnitCell& ucell = env.ucell();
    const int npol = env.npol();
    const int nlocal = pv.get_global_row_size();
    double* const dsloc[3] = {env.fsr().DSloc_x, env.fsr().DSloc_y, env.fsr().DSloc_z};

    const char transN = 'N';
    const char transT = 'T';
    const double one = 1.0;
    const double zero = 0.0;
    assert(nlocal > 0);

    std::vector<double> dm_pot_onsite_dSm(pv.nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        double* tmp_ptr = dsloc[dim];

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

        accumulate_diag_force(pv, ucell, dm_pot_onsite_dSm.data(), dim, force_dftu);

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

        accumulate_onsite_force(env.dftu(), pv, ucell, npol,
                                dm_pot_onsite_dSm.data(), dim, force_dftu);

    } // end dim
    ModuleBase::timer::end("DFTU_LCAO", "cal_force_gamma");
}

void cal_stress_gamma(const DftuFsEnv& env,
                      const double* rho_pot_onsite,
                      ModuleBase::matrix& stress_dftu)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_stress_gamma");
    ModuleBase::timer::start("DFTU_LCAO", "cal_stress_gamma");

    const Parallel_Orbitals& pv = env.pv();
    const UnitCell& ucell = env.ucell();
    const Grid_Driver& gd = env.gd();
    ForceStressArrays& fsr = env.fsr();
    const int npol = env.npol();
    const std::string& ks_solver = env.ks_solver();
    const std::vector<double>& orb_cutoff = env.orb_cutoff();
    const int nlocal = pv.get_global_row_size();
    double* dsloc_x = fsr.DSloc_x;
    double* dsloc_y = fsr.DSloc_y;
    double* dsloc_z = fsr.DSloc_z;
    double* dh_r = fsr.DH_r;

    // shared folding context: read-only params bundled for fold_dSR_gamma
    DFTU_LCAO::FoldingCtx fold_ctx{npol, ks_solver, orb_cutoff, &ucell, &pv, &gd};

    const char transN = 'N';
    const double zero = 0.0;
    const double minus_half = -0.5;

    std::vector<double> dSR_gamma(pv.nloc);
    std::vector<double> dm_pot_onsite_sover(pv.nloc);

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            DFTU_LCAO::fold_dSR_gamma(fold_ctx, dsloc_x, dsloc_y, dsloc_z, dh_r, dim1, dim2, &dSR_gamma[0]);

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

            accumulate_diag_stress(pv, dm_pot_onsite_sover.data(), dim1, dim2, 2.0, stress_dftu);

        } // end dim2
    }     // end dim1
    ModuleBase::timer::end("DFTU_LCAO", "cal_stress_gamma");
}

namespace
{

/// @brief Abort if the caller left required folded-matrix buffers unallocated.
///
/// Force needs the three derivative-overlap arrays; stress additionally
/// needs the derivative-Hamiltonian array. The gamma and multik paths use
/// different buffer sets but identical logic.
void check_folded_arrays(const ForceStressArrays& fsr,
                         const bool cal_force,
                         const bool cal_stress,
                         const bool gamma_only_local)
{
    const double* ds0 = gamma_only_local ? fsr.DSloc_x : fsr.DSloc_Rx;
    const double* ds1 = gamma_only_local ? fsr.DSloc_y : fsr.DSloc_Ry;
    const double* ds2 = gamma_only_local ? fsr.DSloc_z : fsr.DSloc_Rz;
    const bool missing_ds = ds0 == nullptr || ds1 == nullptr || ds2 == nullptr;

    if (cal_force && missing_ds)
    {
        const char* message = gamma_only_local
            ? "fsr.DSloc_x/y/z are nullptr in gamma_only path; the caller must allocate and fill them. "
              "See notes in source/source_lcao/force_stress_lcao.cpp."
            : "fsr.DSloc_Rx/Ry/Rz are nullptr in multik path; the caller must allocate and fill them. "
              "See notes in source/source_lcao/force_stress_lcao.cpp.";
        ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress", message);
    }
    if (cal_stress && (missing_ds || fsr.DH_r == nullptr))
    {
        const char* message = gamma_only_local
            ? "fsr.DSloc_x/y/z or fsr.DH_r is nullptr in gamma_only path; "
              "the caller must allocate and fill them. "
              "See notes in source/source_lcao/force_stress_lcao.cpp."
            : "fsr.DSloc_Rx/Ry/Rz or fsr.DH_r is nullptr in multik path; "
              "the caller must allocate and fill them. "
              "See notes in source/source_lcao/force_stress_lcao.cpp.";
        ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress", message);
    }
}

/// @brief Validate caller-provided folded matrices and the ks_solver layout.
///
/// Fail early with a clear message instead of letting pdgemm_ dereference a
/// nullptr. The folded buffers are column-major local ScaLAPACK blocks, so
/// only column-major ks_solvers are accepted.
void validate_fs_inputs(const DftuFsEnv& env,
                        const bool cal_force,
                        const bool cal_stress,
                        const bool gamma_only_local)
{
    check_folded_arrays(env.fsr(), cal_force, cal_stress, gamma_only_local);

    if ((cal_force || cal_stress)
        && !ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(env.ks_solver()))
    {
        ModuleBase::WARNING_QUIT("DFTU_LCAO::force_stress",
            "non column-major ks_solver is not supported for DFT+U force/stress; "
            "the folded matrix layout assumption would be violated");
    }
}

/// @brief Build rho * V_onsite and accumulate gamma-point force/stress over k points.
void run_gamma_loop(const DftuFsEnv& env,
                    const bool cal_force,
                    const bool cal_stress,
                    const std::vector<std::vector<double>>* dmk_d,
                    const K_Vectors& kv,
                    ModuleBase::matrix& force_dftu,
                    ModuleBase::matrix& stress_dftu)
{
    assert(dmk_d != nullptr);
    Plus_U_Base& dftu = env.dftu();
    const UnitCell& ucell = env.ucell();
    const Parallel_Orbitals& pv = env.pv();
    const int npol = env.npol();
    const int nlocal = pv.get_global_row_size();

    const char transN = 'N';
    const char transT = 'T';
    const double alpha = 1.0;
    const double beta = 0.0;

    std::vector<double> rho_pot_onsite(pv.nloc);

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        const int spin = kv.isk[ik];
        std::vector<double> pot_onsite(pv.nloc, 0.0);

        pot_onsite_real(dftu, ucell, &pv, spin, false, pot_onsite.data(), npol);

#ifdef __MPI
        ScalapackConnector::gemm(transT, transN, nlocal, nlocal, nlocal,
                alpha, (*dmk_d)[spin].data(), 1, 1,
                pv.desc, pot_onsite.data(), 1, 1,
                pv.desc, beta, rho_pot_onsite.data(),
                1, 1, pv.desc);
#endif

        if (cal_force)
        {
            cal_force_gamma(env, rho_pot_onsite.data(), force_dftu);
        }
        if (cal_stress)
        {
            cal_stress_gamma(env, rho_pot_onsite.data(), stress_dftu);
        }
    } // ik
}

/// @brief Build rho * V_onsite and accumulate multik force/stress over k points.
void run_k_loop(const DftuFsEnv& env,
                const bool cal_force,
                const bool cal_stress,
                const std::vector<std::vector<std::complex<double>>>* dmk_c,
                const K_Vectors& kv,
                ModuleBase::matrix& force_dftu,
                ModuleBase::matrix& stress_dftu)
{
    assert(dmk_c != nullptr);
    Plus_U_Base& dftu = env.dftu();
    const UnitCell& ucell = env.ucell();
    const Parallel_Orbitals& pv = env.pv();
    const int npol = env.npol();
    const int nlocal = pv.get_global_row_size();

    const char transN = 'N';
    const char transT = 'T';
    const int one_int = 1;
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double> beta(0.0, 0.0);

    std::vector<std::complex<double>> rho_pot_onsite(pv.nloc);

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        const int spin = kv.isk[ik];
        std::vector<std::complex<double>> pot_onsite(pv.nloc, std::complex<double>(0.0, 0.0));

        pot_onsite_complex(dftu, ucell, &pv, spin, false, pot_onsite.data(), npol);

#ifdef __MPI
        ScalapackConnector::gemm(transT, transN, nlocal, nlocal, nlocal,
                alpha, (*dmk_c)[ik].data(), one_int, one_int,
                pv.desc, pot_onsite.data(), one_int, one_int, pv.desc, beta,
                rho_pot_onsite.data(), one_int, one_int, pv.desc);
#endif

        if (cal_force)
        {
            cal_force_k(env, ik, kv.kvec_d[ik], rho_pot_onsite.data(), force_dftu);
        }
        if (cal_stress)
        {
            cal_stress_k(env, ik, kv.kvec_d[ik], rho_pot_onsite.data(), stress_dftu);
        }
    } // ik
}

/// @brief MPI-reduce the stress, mirror the upper triangle and apply cell scaling.
void finalize_stress(const UnitCell& ucell, ModuleBase::matrix& stress_dftu)
{
    Parallel_Reduce::reduce_pool(stress_dftu.c, stress_dftu.nr * stress_dftu.nc);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i > j)
            {
                stress_dftu(i, j) = stress_dftu(j, i);
            }
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

} // namespace

void force_stress(const DftuFsEnv& env,
                  const bool cal_force,
                  const bool cal_stress,
                  std::vector<std::vector<double>>* dmk_d,
                  std::vector<std::vector<std::complex<double>>>* dmk_c,
                  ModuleBase::matrix& force_dftu,
                  ModuleBase::matrix& stress_dftu,
                  const K_Vectors& kv,
                  const bool gamma_only_local)
{
    ModuleBase::TITLE("DFTU_LCAO", "force_stress");
    ModuleBase::timer::start("DFTU_LCAO", "force_stress");

    validate_fs_inputs(env, cal_force, cal_stress, gamma_only_local);

    if (cal_force)
    {
        force_dftu.zero_out();
    }
    if (cal_stress)
    {
        stress_dftu.zero_out();
    }

    if (gamma_only_local)
    {
        run_gamma_loop(env, cal_force, cal_stress, dmk_d, kv, force_dftu, stress_dftu);
    }
    else
    {
        run_k_loop(env, cal_force, cal_stress, dmk_c, kv, force_dftu, stress_dftu);
    }

    if (cal_force)
    {
        Parallel_Reduce::reduce_pool(force_dftu.c, force_dftu.nr * force_dftu.nc);
    }
    if (cal_stress)
    {
        finalize_stress(env.ucell(), stress_dftu);
    }

    ModuleBase::timer::end("DFTU_LCAO", "force_stress");
}

} // namespace DFTU_LCAO
