
namespace LCAO_domain
{

void init_basis_lcao(const Parallel_Orbitals& pv,
		const Input_para& inp, 
		UnitCell& ucell,
        TwoCenterBundle& two_center_bundle)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "init_basis_lcao");

    const int nlocal = GlobalV::NLOCAL;

    // autoset NB2D first
    if (GlobalV::NB2D == 0)
    {
        if (nlocal > 0)
        {
            GlobalV::NB2D = (GlobalV::NSPIN == 4) ? 2 : 1;
        }
        if (nlocal > 500)
        {
            GlobalV::NB2D = 32;
        }
        if (nlocal > 1000)
        {
            GlobalV::NB2D = 64;
        }
    }

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.

    two_center_bundle.build_orb(ucell.ntype, ucell.orbital_fn);
    two_center_bundle.build_alpha(GlobalV::deepks_setorb, &ucell.descriptor_file);
    two_center_bundle.build_orb_onsite(PARAM.inp.onsite_radius);
    // currently deepks only use one descriptor file, so cast bool to int is
    // fine

    // TODO Due to the omnipresence of GlobalC::ORB, we still have to rely
    // on the old interface for now.
    two_center_bundle.to_LCAO_Orbitals(GlobalC::ORB, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);

    ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, GlobalC::ORB);

    two_center_bundle.build_beta(ucell.ntype, ucell.infoNL.Beta);

    int Lmax = 0;
#ifdef __EXX
    Lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

#ifdef USE_NEW_TWO_CENTER
    two_center_bundle.tabulate();
#else
    two_center_bundle.tabulate(inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax);
#endif

    // setup_2d_division
#ifdef __MPI
    // storage form of H and S matrices on each processor
    // is determined in 'divide_HS_2d' subroutine

    int try_nb = pv.init(nlocal, nlocal, GlobalV::NB2D, DIAG_WORLD);
    try_nb += pv.set_nloc_wfc_Eij(GlobalV::NBANDS, GlobalV::ofs_running, GlobalV::ofs_warning);
    if (try_nb != 0)
    {
        pv.set(nlocal, nlocal, 1, pv.blacs_ctxt);
        try_nb = pv.set_nloc_wfc_Eij(GlobalV::NBANDS, GlobalV::ofs_running, GlobalV::ofs_warning);
    }

    // init blacs context for genelpa
    pv.set_desc_wfc_Eij(nlocal, GlobalV::NBANDS, pv.nrow);

#else
    pv.set_serial(nlocal, nlocal);
    pv.nrow_bands = nlocal;
    pv.ncol_bands = GlobalV::NBANDS;
    // Zhang Xiaoyang enable the serial version of LCAO and recovered this function usage. 2024-07-06
#endif

    pv.set_atomic_trace(GlobalC::ucell.get_iat2iwt(), GlobalC::ucell.nat, nlocal);

    return;
}

}
