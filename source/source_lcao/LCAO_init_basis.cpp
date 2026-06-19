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

    two_center_bundle.build_orb(ucell.ntype, ucell.orbital_fn.data(), PARAM.inp.orbital_dir);
    two_center_bundle.build_alpha(PARAM.globalv.deepks_setorb, &ucell.descriptor_file);
    two_center_bundle.build_orb_onsite(onsite_radius);
    // currently deepks only use one descriptor file, so cast bool to int is
    // fine

    // TODO Due to the omnipresence of LCAO_Orbitals, we still have to rely
    // on the old interface for now.
    two_center_bundle.to_LCAO_Orbitals(orb, lcao_ecut, lcao_dk, lcao_dr, lcao_rmax,
                                       PARAM.inp.out_element_info, PARAM.inp.cal_force);

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
    // (p <= q) with square block nb2d. The time-vs-nb2d curve is U-shaped: too large
    // -> load-imbalance cliff (one process owns a whole panel); too small -> poor BLAS
    // and heavy block-cyclic communication. Healthy window [nb_lo, nb_hi] below.
    // Only scalapack_gvx diagonalizes on this 2D grid (genelpa/elpa tune internally;
    // lapack/cusolver/pexsi do not), so the check is restricted to it.
    if (PARAM.inp.ks_solver == "scalapack_gvx")
    {
        const int kpar = (PARAM.globalv.kpar_lcao > 0) ? PARAM.globalv.kpar_lcao : 1;
        // Processes running one (per-pool) diagonalization:
        //   kpar == 1 : the grid is this ParaV grid (built on DIAG_WORLD) -> pv.dim0*pv.dim1.
        //   kpar  > 1 : hsolver re-splits MPI_COMM_WORLD into kpar pools of NPROC/kpar ranks
        //               (not the DIAG_WORLD grid). Uneven pools (NPROC % kpar != 0) are skipped.
        int np_pool = 0; // processes per pool (0 => skip the check)
        if (kpar <= 1)
        {
            np_pool = pv.dim0 * pv.dim1;
        }
        else if (GlobalV::NPROC % kpar == 0)
        {
            np_pool = GlobalV::NPROC / kpar;
        }
        if (np_pool > 1 && nlocal > 0)
        {
            // near-square factorization np_pool = p * q, p <= q (matches Parallel_2D)
            int p_row = static_cast<int>(std::sqrt(np_pool + 0.5));
            while (p_row > 1 && np_pool % p_row != 0) { --p_row; }
            const int p_col = np_pool / p_row;    // long edge q (>= p_row)

            // Two-sided window: nb_hi = floor(N/2q) keeps >= 2 blocks per process;
            // nb_lo = min(16, nb_hi) avoids tiny blocks; recommended = min(64, nb_hi).
            // nspin==4 carries 2-component spinors that must stay paired in one block
            // (hence autoset/fallback use nb2d=2, not 1), so snap the window to a multiple
            // of 2 -- an odd nb2d would break the spinor blocking and segfault.
            const int nb_unit = (PARAM.inp.nspin == 4) ? 2 : 1;
            auto snap = [nb_unit](int v) { v = v / nb_unit * nb_unit; return v < nb_unit ? nb_unit : v; };
            const int nb_hi = snap((nlocal >= 2 * p_col) ? nlocal / (2 * p_col) : 1);
            const int nb_lo = snap((16 < nb_hi) ? 16 : nb_hi);
            const int nb_opt = snap((nb_hi < 64) ? nb_hi : 64);
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

            // user-set nb2d (!=0): keep the value, only warn. auto nb2d (==0): correct it
            // to nb_opt (pv.nb feeds both the kpar==1 path and the per-pool Parallel_K2D).
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
                    // Validate nb_opt like the initial distribution: set_nloc_wfc_Eij
                    // returns non-zero if it is incompatible with the band/grid layout
                    // (ceil(nbands/nb_opt) < grid width). If so, revert to the validated
                    // nb_cur -- a half-updated pv would crash the later wavefunction setup.
                    int retry = pv.set(nlocal, nlocal, nb_opt, pv.blacs_ctxt);
                    retry += pv.set_nloc_wfc_Eij(PARAM.inp.nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
                    if (retry != 0)
                    {
                        pv.set(nlocal, nlocal, nb_cur, pv.blacs_ctxt);
                        pv.set_nloc_wfc_Eij(PARAM.inp.nbands, GlobalV::ofs_running, GlobalV::ofs_warning);
                        GlobalV::ofs_warning << "init_basis_lcao: auto nb2d=" << nb_cur << " is " << issue
                            << " for N=" << nlocal << ", kpar=" << kpar << "; recommended nb2d=" << nb_opt
                            << " is incompatible with the band/grid layout, so nb2d=" << nb_cur << " is kept.\n";
                    }
                    else
                    {
                        GlobalV::ofs_warning << "init_basis_lcao: auto nb2d=" << nb_cur << " is " << issue
                            << " for N=" << nlocal << ", kpar=" << kpar << "; auto-adjusted to nb2d=" << nb_opt << ".\n";
                    }
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
