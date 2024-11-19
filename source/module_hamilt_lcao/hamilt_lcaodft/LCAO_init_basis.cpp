#include "LCAO_domain.h"

#include "module_parameter/parameter.h"
/// once the GlobalC::exx_info has been deleted, this include can be gone 
/// mohan note 2024-07-21
#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif

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

    two_center_bundle.build_orb(ucell.ntype, ucell.orbital_fn);
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

    int Lmax = 0;
#ifdef __EXX
    Lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

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
