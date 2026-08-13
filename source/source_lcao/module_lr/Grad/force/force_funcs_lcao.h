#pragma once
#include "source_pw/module_pwdft/force_pw.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_io/module_output/output_log.h"

// This file extracts out some terms from the ground state force code 
// to calculate the excited state force with a different density matrix.
template<typename T>
class ForcePWTerms
{
public:
    ModuleBase::matrix operator()(const UnitCell& ucell,
        const Charge& chr,
        const ModulePW::PW_Basis& rhopw,
        const pseudopot_cell_vl& locpp,
        const Structure_Factor& sf,
        const bool with_ewald = true,
        const elecstate::ElecState* pelec = nullptr)
    {
        if (PARAM.inp.nspin == 4) { throw std::runtime_error("ForcePWTerms: nspin=4 is not supported."); }
        ModuleBase::TITLE("Force_Stress_LCAO", "cal_force_pw");
        Forces<T> f_pw(ucell.nat);
        ModuleBase::matrix fvl_dvl(ucell.nat, 3), fewalds(ucell.nat, 3), fcc(ucell.nat, 3), fscc(ucell.nat, 3);
        //--------------------------------------------------------
        // local pseudopotential force:
        // use charge density; plane wave; local pseudopotential;
        //--------------------------------------------------------
        f_pw.cal_force_loc(ucell, fvl_dvl, &rhopw, locpp.vloc, &chr);
        //--------------------------------------------------------
        // ewald force: use plane wave only.
        //--------------------------------------------------------
        if (with_ewald)
        {
            f_pw.cal_force_ew(ucell, fewalds, &rhopw, &sf);
        }
        //--------------------------------------------------------
        // force due to core correlation.
        //--------------------------------------------------------
        UnitCell& ucell_noconst = const_cast<UnitCell&>(ucell);
        f_pw.cal_force_cc(fcc, &rhopw, &chr, locpp.numeric, ucell_noconst); // no problem for nspin=1 and 2
        //--------------------------------------------------------
        // force due to self-consistent charge (invalid in from-scratch LR case)
        //--------------------------------------------------------
        if (pelec)
        {
            f_pw.cal_force_scc(fscc, &rhopw, pelec->vnew, pelec->vnew_exist, locpp.numeric, ucell);
        }
        if (PARAM.inp.test_force)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "VL_dVL      FORCE (eV/Angstrom)", fvl_dvl, false);
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "EWALD      FORCE (eV/Angstrom)", fewalds, false);
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "NLCC      FORCE (eV/Angstrom)", fcc, false);
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "SCC      FORCE (eV/Angstrom)", fscc, false);
        }
        return fvl_dvl + fewalds + fcc + fscc;
    }
};

// calculate force and stress for Nonlocal part
// nspin = 2 or 4 will not be used in the current excited state force calculation
// for nspin = 1 or 2
template<typename TK>
ModuleBase::matrix cal_force_nonlocal(
    const UnitCell& ucell,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    const Grid_Driver& gd,
    const TwoCenterBundle& two_center_bundle,
    const elecstate::DensityMatrix<TK, double>& dm    )
{
    ModuleBase::TITLE("Force_Stress_LCAO", "cal_force_nonlocal_dvnl");
    std::vector<double> orb_cutoffs(ucell.ntype);
    for (int it = 0;it < ucell.ntype;++it) { orb_cutoffs[it] = ucell.atoms[it].Rcut; }

    hamilt::Nonlocal<hamilt::OperatorLCAO<TK, double>> tmp_nonlocal(nullptr,
        kvec_d,
        nullptr,
        &ucell,
        orb_cutoffs,
        &gd,
        two_center_bundle.overlap_orb_beta.get());
    const int nspin = dm.get_DMR_vector().size();
    if(nspin==2)
    {
        const_cast<elecstate::DensityMatrix<TK, double>*>(&dm)->switch_dmr(1); //spin-up + spin-down
    }
    const hamilt::HContainer<double>* dmr = dm.get_DMR_pointer(1);
    ModuleBase::matrix fvnl(ucell.nat, 3);
    ModuleBase::matrix svnl; // no use now, only for passing into interfaces
    tmp_nonlocal.cal_force_stress(/*force*/true, /*stress*/false, dmr, fvnl, svnl);
    if (nspin == 2)
    {
        const_cast<elecstate::DensityMatrix<TK, double>*>(&dm)->switch_dmr(0);
    }
    return fvnl;
}
// calculate force and stress for Nonlocal part (Hellmann-Feynman term)
// for nspin = 4
template<typename TK>
ModuleBase::matrix cal_force_nonlocal_dvnl(
    UnitCell ucell,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    const Grid_Driver& gd,
    const TwoCenterBundle& two_center_bundle,
    const elecstate::DensityMatrix<TK, std::complex<double>>& dm)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "cal_force_nonlocal_dvnl");

    std::vector<double> orb_cutoffs(ucell.ntype);
    for (int it = 0;it < ucell.ntype;++it) { orb_cutoffs[it] = ucell.atoms[it].Rcut; }

    hamilt::Nonlocal<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_nonlocal(
        nullptr,
        kvec_d,
        nullptr,
        &ucell,
        orb_cutoffs,
        &gd,
        two_center_bundle.overlap_orb_beta.get());

    hamilt::HContainer<std::complex<double>> tmp_dmr(dm.get_DMR_pointer(1)->get_paraV());
    std::vector<int> ijrs = dm.get_DMR_pointer(1)->get_ijr_info();
    tmp_dmr.insert_ijrs(&ijrs);
    tmp_dmr.allocate();
    dm.cal_DMR_full(&tmp_dmr);
    ModuleBase::matrix fvnl(ucell.nat, 3);
    ModuleBase::matrix svnl; // no use now, only for passing into interfaces
    tmp_nonlocal.cal_force_stress(/*force*/true, /*stress*/false, &tmp_dmr, fvnl, svnl);
    return fvnl;
}