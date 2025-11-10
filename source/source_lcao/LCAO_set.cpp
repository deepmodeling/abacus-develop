#include "source_lcao/LCAO_set_pot.h"
#include "source_io/module_parameter/parameter.h"

namespace LCAO_domain
{

void set_psi_occ_dm_chg(psi, kv, pv, inp, pelec, dmat, chr)
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
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "read electronic wave functions failed");
        }
    }

    // 3) set occupations, tddft does not need to set occupations in the first scf
    if (inp.ocp && inp.esolver_type != "tddft")
    {
        elecstate::fixed_weights(inp.ocp_kb, inp.nbands, inp.nelec,
          pelec->klist, pelec->wg, pelec->skip_weights);
    }

    // 4) init DMK, but DMR is constructed in before_scf()
    dmat.allocate_dm(&kv, &pv, inp.nspin);

    // 5) init charge density
    chr.allocate(inp.nspin);
}


void set_pot(pelec, locpp, pw_rho, pw_rhod, sf, solvent, dftu, pv, kv, orb_, exx_nao, deepks)
{
    // 1) init local pseudopotentials
    locpp.init_vloc(ucell, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    // 2) init potentials
    if (pelec->pot == nullptr)
    {
        // where is the pot deleted?
        pelec->pot = new elecstate::Potential(pw_rhod, pw_rho,
          &ucell, &(locpp.vloc), &(sf), &(solvent),
          &(pelec->f_en.etxc), &(pelec->f_en.vtxc));
    }

    // 3) initialize DFT+U
    if (inp.dft_plus_u)
    {
        dftu.init(ucell, &pv, kv.get_nks(), &orb_);
    }

    // 4) init exact exchange calculations
    exx_nao.before_runner(ucell, kv, orb_, pv, inp);

    // 5) init deepks
    deepks.before_runner(ucell, kv.get_nks(), orb_, pv, inp);

    return;
}


}
