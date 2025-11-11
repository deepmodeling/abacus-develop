#include "source_lcao/LCAO_set_pot.h"
#include "source_io/module_parameter/parameter.h"

template <typename TK>
void LCAO_domain::set_psi_occ_dm_chg(
		const K_Vectors &kv, // k-points
		psi::Psi<TK>* psi, // coefficients of NAO basis
		const Parallel_Orbitals &pv, // parallel scheme of NAO basis
		elecstate::ElecState* pelec, // eigen values and weights
		LCAO_domain::Setup_DM<TK> &dmat, // density matrix 
		Charge &chr, // charge density 
		const Input_para &inp) // input parameters
{
    //! 1) init electronic wave function psi
    Setup_Psi<TK>::allocate_psi(psi, kv, pv, inp);

    //! 2) read psi from file
    if (inp.init_wfc == "file" && inp.esolver_type != "tddft")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
             pv, *psi, pelec->ekb, pelec->wg, kv.ik2iktot,
             kv.get_nkstot(), inp.nspin))
        {
            ModuleBase::WARNING_QUIT("set_psi_occ_dm_chg", "read electronic wave functions failed");
        }
    }

    //! 3) set occupations, tddft does not need to set occupations in the first scf
    if (inp.ocp && inp.esolver_type != "tddft")
    {
        elecstate::fixed_weights(inp.ocp_kb, inp.nbands, inp.nelec,
          kv, pelec->wg, pelec->skip_weights);
    }

    //! 4) init DMK, but DMR is constructed in before_scf()
    dmat.allocate_dm(&kv, &pv, inp.nspin);

    //! 5) init charge density
    chr.allocate(inp.nspin);

    return;
}



void LCAO_domain::set_pot(
		const K_Vectors &kv, 
	    const Structure_Factor& sf,	
		const ModulePW::PW_Basis &pw_rho, 
		const ModulePW::PW_Basis &pw_rhod, 
		elecstate::ElecState* pelec,
		const LCAO_Orbitals& orb,
		const Parallel_Orbitals &pv, 
		pseudopot_cell_vl &locpp, 
        Plus_U &dftu,
        surchem& solvent,
        Exx_NAO<T> &exx_nao,
        Setup_DeePKS<T> &deepks,
        const Input_para &inp)
{
    //! 1) init local pseudopotentials
    locpp.init_vloc(ucell, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //! 2) init potentials
    if (pelec->pot == nullptr)
    {
        // where is the pot deleted?
        pelec->pot = new elecstate::Potential(&pw_rhod, &pw_rho,
          &ucell, &(locpp.vloc), &(sf), &(solvent),
          &(pelec->f_en.etxc), &(pelec->f_en.vtxc));
    }

    //! 3) initialize DFT+U
    if (inp.dft_plus_u)
    {
        dftu.init(ucell, &pv, kv.get_nks(), &orb);
    }

    //! 4) init exact exchange calculations
    exx_nao.before_runner(ucell, kv, orb, pv, inp);

    //! 5) init deepks
    deepks.before_runner(ucell, kv.get_nks(), orb, pv, inp);

    return;
}


