#include "LCAO_domain.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/parallel_comm.h"

#include <cmath>

namespace LCAO_domain
{

void init_basis_lcao(Parallel_Orbitals& pv,
        const double &onsite_radius,
        const double &lcao_ecut,
        const double &lcao_dk,
        const double &lcao_dr,
        const double &lcao_rmax,
		UnitCell& ucell,
        TwoCenterBundle& two_center_bundle,
        LCAO_Orbitals& orb
)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_basis_lcao");

    const int nlocal = PARAM.globalv.nlocal;
    int nb2d = PARAM.inp.nb2d;
    // autoset NB2D first
    if (nb2d == 0)
    {
        if (nlocal > 0)
        {
            nb2d = (PARAM.inp.nspin == 4) ? 2 : 1;
        }
        if (nlocal > 500)
        {
            nb2d = 32;
        }
        if (nlocal > 1000)
        {
            nb2d = 64;
        }
    }

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.

    two_center_bundle.build_orb(ucell.ntype, ucell.orbital_fn.data());
    two_center_bundle.build_alpha(PARAM.globalv.deepks_setorb, &ucell.descriptor_file);
    two_center_bundle.build_orb_onsite(onsite_radius);
    // currently deepks only use one descriptor file, so cast bool to int is
    // fine

    // TODO Due to the omnipresence of LCAO_Orbitals, we still have to rely
    // on the old interface for now.
    two_center_bundle.to_LCAO_Orbitals(orb, lcao_ecut, lcao_dk, lcao_dr, lcao_rmax);

    if (PARAM.inp.vnl_in_h)
    {
        ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, orb);
        two_center_bundle.build_beta(ucell.ntype, ucell.infoNL.Beta);
    }

#ifdef USE_NEW_TWO_CENTER
    two_center_bundle.tabulate();
#else
    two_center_bundle.tabulate(lcao_ecut, lcao_dk, lcao_dr, lcao_rmax);
#endif

    // setup_2d_division
#ifdef __MPI
    // storage form of H and S matrices on each processor
    // is determined in 'divide_HS_2d' subroutine

    int try_nb = pv.init(nlocal, nlocal, nb2d, DIAG_WORLD);
    try_nb += pv.set_nloc_wfc_Eij(PARAM.inp.nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
    if (try_nb != 0)
    {
        // fall back to the minimum size, 1 or 2 (nspin=4)
        const int min_size = (PARAM.inp.nspin == 4) ? 2 : 1;
        pv.set(nlocal, nlocal, min_size, pv.blacs_ctxt);
        try_nb = pv.set_nloc_wfc_Eij(PARAM.inp.nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
    }

    // ---- nb2d (ScaLAPACK 2D block-cyclic block size) load-balance check ----
    // ScaLAPACK diagonalizes the N x N matrix (N = nlocal) on a p x q process grid
    // (p <= q) with square block size nb2d. With kpar k-point pools the diagonalization
    // runs *per pool* on NPROC/kpar processes (Parallel_K2D::P2D_pool, whose block size
    // is ParaV->get_block_size() == pv.nb), so the effective grid is the near-square
    // factorization of NPROC/kpar -- exactly how Parallel_2D builds it. The long edge q
    // governs load balance:
    //   B = N / (nb2d * q)   blocks owned per process along q.
    //   B <= 1 : one process owns a whole panel while others idle -> a catastrophic
    //            load-imbalance "cliff" (nb2d too large; energy unaffected, time blows up).
    // But nb2d also must not be too SMALL: blocks below ~16 lose BLAS/GEMM efficiency and
    // explode block-cyclic communication (nb=1 = the slow "over-scatter" end). So the time-
    // vs-nb2d curve is U-shaped and the healthy window is two-sided:
    //     nb_lo = min(16, N/(2q))  <=  nb2d  <=  nb_hi = floor(N/(2q))
    // recommended nb2d = min(64, nb_hi) (largest balanced block, capped for BLAS).
    // This block-size U-curve is a property of the ScaLAPACK 2D block-cyclic *dense*
    // diagonalization, so the check is restricted to ks_solver=scalapack_gvx. (genelpa/
    // elpa do their own internal block tuning; lapack/cusolver/pexsi do not diagonalize
    // on this distributed 2D grid, so the nb2d cliff does not apply to them.)
    if (PARAM.inp.ks_solver == "scalapack_gvx")
    {
        const int kpar = (PARAM.globalv.kpar_lcao > 0) ? PARAM.globalv.kpar_lcao : 1;
        const int np_total = pv.dim0 * pv.dim1;   // diagonalization world size
        // Pools split the world evenly; if np_total is not divisible by kpar the
        // pool size is ill-defined, so skip the (purely advisory) heuristic rather
        // than estimate the grid from a truncated np_total/kpar.
        const bool pool_well_defined = (kpar > 0 && np_total % kpar == 0);
        const int np_pool = pool_well_defined ? np_total / kpar : 0; // processes per pool
        if (np_pool > 1 && nlocal > 0)
        {
            // near-square factorization np_pool = p * q, p <= q (matches Parallel_2D)
            int p_row = static_cast<int>(std::sqrt(np_pool + 0.5));
            while (p_row > 1 && np_pool % p_row != 0) { --p_row; }
            const int p_col = np_pool / p_row;    // long edge q (>= p_row)

            // Healthy block-size window is two-sided (the nb2d-vs-time curve is U-shaped):
            //   nb_hi = floor(N / (2*q))  upper bound -- keep >= 2 blocks per process
            //                             (B >= 2); larger nb -> load imbalance "cliff".
            //   nb_lo = min(16, nb_hi)    lower bound -- blocks below ~16 lose BLAS/GEMM
            //                             efficiency (1x1 ops) and explode block-cyclic
            //                             communication (panel count ~ N/nb); nb=1 is the
            //                             slow "over-scatter" end. Capped by nb_hi because a
            //                             tiny system on many processes cannot afford large
            //                             blocks (then balance wins and nb_lo == nb_hi).
            // recommended = the largest balanced block, capped at 64 for BLAS efficiency.
            const int nb_hi = (nlocal >= 2 * p_col) ? nlocal / (2 * p_col) : 1;
            const int nb_lo = (16 < nb_hi) ? 16 : nb_hi;
            const int nb_opt = (nb_hi < 64) ? nb_hi : 64;
            const int nb_cur = pv.nb;

            const char* issue = nullptr;
            if (nb_cur > nb_hi)
            {
                issue = "too large -> ScaLAPACK load imbalance (one process owns a whole panel)";
            }
            else if (nb_cur < nb_lo)
            {
                issue = "too small -> over-scatter (poor BLAS efficiency and heavy communication)";
            }

            // Two cases, both reported via ofs_warning:
            //   (1) user explicitly set nb2d (PARAM.inp.nb2d != 0): respect it -- do NOT
            //       change the value -- but warn so the issue is visible.
            //   (2) auto nb2d (PARAM.inp.nb2d == 0): correct it to nb_opt. pv.nb feeds both
            //       the kpar==1 path and Parallel_K2D (ParaV->get_block_size()), so
            //       resetting it here also fixes the per-pool diagonalization.
            if (issue != nullptr)
            {
                if (PARAM.inp.nb2d != 0)
                {
                    GlobalV::ofs_warning << "init_basis_lcao: user-set nb2d=" << nb_cur << " is " << issue
                        << " for N=" << nlocal << ", kpar=" << kpar << " (per-pool grid " << p_row << "x"
                        << p_col << "); recommended nb2d=" << nb_opt << " (user value kept, not changed).\n";
                }
                else
                {
                    GlobalV::ofs_warning << "init_basis_lcao: auto nb2d=" << nb_cur << " is " << issue
                        << " for N=" << nlocal << ", kpar=" << kpar << "; auto-adjusted to nb2d=" << nb_opt << ".\n";
                    pv.set(nlocal, nlocal, nb_opt, pv.blacs_ctxt);
                    pv.set_nloc_wfc_Eij(PARAM.inp.nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
                }
            }
        }
    }

    // init blacs context for genelpa
    pv.set_desc_wfc_Eij(nlocal, PARAM.inp.nbands, pv.nrow);

#else
    pv.set_serial(nlocal, nlocal);
    pv.nrow_bands = nlocal;
    pv.ncol_bands = PARAM.inp.nbands;
    // Zhang Xiaoyang enable the serial version of LCAO and recovered this function usage. 2024-07-06
#endif

    pv.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, nlocal);

    return;
}

}
