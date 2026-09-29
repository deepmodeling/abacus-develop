#include "force_stress_lcao.h"

#include "force_stress_assemble.h"
#include "force_stress_terms.h"

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_base/parallel_reduce.h"
#include "source_pw/module_pwdft/dftu_base.h" //Quxin add for DFT+U on 20201029
#include "source_lcao/module_dftu/dftu_nao_fs_k.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
// new
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_pot/h_tddft_pw.h"       // Taoni add 2025-02-20
#include "source_estate/module_pot/efield.h"           // liuyu add 2022-05-18
#include "source_estate/module_pot/gatefield.h"        // liuyu add 2022-09-13
#include "source_hamilt/module_surchem/surchem.h" //sunml add 2022-08-10
#include "source_hamilt/module_vdw/vdw.h"
#ifdef __MLALGO
#include "source_lcao/module_deepks/lcao_deepks.h"    //caoyu add for deepks 2021-06-03
#include "source_lcao/module_deepks/lcao_deepks_io.h" // mohan add 2024-07-22
#include "source_lcao/module_deepks/deepks_force.h"
#endif
#include "source_lcao/module_dftu/dftu_nao_adj.h"
#include "source_lcao/module_dftu/dftu_nao_fs_r.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/overlap.h"
#include "source_lcao/module_operator_lcao/td_pot_hybrid.h"
#include "source_lcao/pulay_fs.h"
#include "source_lcao/module_rt/force_rt_overlap.h"


// mohan add 2025-11-04
template <>
void assign_dmk_ptr<double>(
    module_dm::DensityMatrix<double,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c
) {
    std::vector<std::vector<double>>& dmk_tmp = dm->get_dmk_vec();
    dmk_d = &dmk_tmp;
    dmk_c = nullptr;
}

template <>
void assign_dmk_ptr<std::complex<double>>(
    module_dm::DensityMatrix<std::complex<double>,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c
) {
    std::vector<std::vector<std::complex<double>>>& dmk_tmp = dm->get_dmk_vec();
    dmk_c = &dmk_tmp;
    dmk_d = nullptr;
}



template <typename T>
Force_Stress_LCAO<T>::Force_Stress_LCAO(Record_adj& ra, const int nat_in) : RA(&ra), nat(nat_in)
{
}
template <typename T>
Force_Stress_LCAO<T>::~Force_Stress_LCAO()
{
}
template <typename T>
void Force_Stress_LCAO<T>::getForceStress(UnitCell& ucell,
                                          const vdw::VdwResult* vdw_result,
                                          const bool isforce,
                                          const bool isstress,
                                          const bool istestf,
                                          const bool istests,
                                          const Grid_Driver& gd,
                                          Parallel_Orbitals& pv,
                                          const elecstate::ElecState* pelec,
                                          module_dm::Setup_DM<T> &dmat, // mohan add 2025-11-03
                                          const psi::Psi<T>* psi,
                                          const TwoCenterBundle& two_center_bundle,
                                          const LCAO_Orbitals& orb,
                                          ModuleBase::matrix& fcs,
                                          ModuleBase::matrix& scs,
                                          const pseudopot_cell_vl& locpp,
                                          const Structure_Factor& sf,
                                          const K_Vectors& kv,
                                          ModulePW::PW_Basis* rhopw,
                                          surchem& solvent,
                                          Plus_U_Base &dftu, // mohan add 2025-11-07
                                          Setup_DeePKS<T>& deepks,
                                          Exx_NAO<T> &exx_nao,
                                          ModuleSymmetry::Symmetry* symm,
                                          const Exx_Info& exx_info,
                                          const FSCalcConfig& cfg,
                                          const int td_stype,
                                          hamilt::Hamilt<T>* p_hamilt)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "getForceStress");
    ModuleBase::timer::start("Force_Stress_LCAO", "getForceStress");

    if (!isforce && !isstress)
    {
        ModuleBase::timer::end("Force_Stress_LCAO", "getForceStress");
        return;
    }

    const int nat = ucell.nat;

    // NOTE: ForceStressArrays is no longer needed as we use operator-based force calculation
    // ForceStressArrays fsr; // removed - no longer needed

    // total force : ModuleBase::matrix fcs;

    // part of total force / stress, grouped so the assembly helpers can take a
    // single container reference.
    LCAOForceParts parts;
    LCAOStressParts sparts;

    parts.fvl_dphi.create(nat, 3); // must do it now, update it later, noted by zhengdy

    if (isforce)
    {
        fcs.create(nat, 3);
        parts.foverlap.create(nat, 3); // overlap force
        parts.ftvnl_dphi.create(nat, 3); // pulay force of NAO
        parts.fvnl_dbeta.create(nat, 3); // pulay force of non-local projectors
        parts.fvl_dvl.create(nat, 3); // force from local potentials
        parts.fewalds.create(nat, 3); // Ewald force
        parts.fcc.create(nat, 3); // force due to core correction
        parts.fscc.create(nat, 3); // force due to self-consistent field
        parts.fvnl_dalpha.create(nat, 3); // deepks
        parts.fpothybrid.create(nat, 3); // pulay force for hybrid gauge rt-tddft

        // calculate basic terms in Force, same method with PW base
        this->calForcePwPart(ucell, parts.fvl_dvl, parts.fewalds, parts.fcc, parts.fscc,
                             pelec->f_en.etxc, pelec->vnew, pelec->vnew_exist, pelec->charge, rhopw,
                             locpp, sf, cfg);
    }

    // total stress : ModuleBase::matrix scs

    //! stress
    if (isstress)
    {
        scs.create(3, 3);
        sparts.sigmacc.create(3, 3);
        sparts.sigmadvl.create(3, 3);
        sparts.sigmaewa.create(3, 3);
        sparts.sigmaxc.create(3, 3);
        sparts.sigmahar.create(3, 3);

        sparts.soverlap.create(3, 3);
        sparts.stvnl_dphi.create(3, 3);
        sparts.svnl_dbeta.create(3, 3);
        sparts.svl_dphi.create(3, 3);
        sparts.svnl_dalpha.create(3, 3);

        // calculate basic terms in Stress, similar method with PW base
        this->sc_pw.stress_pw_terms(ucell, sparts.sigmadvl, sparts.sigmahar, sparts.sigmaewa,
                                    sparts.sigmacc, sparts.sigmaxc, pelec->f_en.etxc, pelec->charge,
                                    rhopw, locpp, sf, cfg.nspin, cfg.domag, cfg.domag_z, cfg.gga_grad, cfg.gamma_only_pw);
    }
    // Calculate operator-based force/stress terms (kinetic, overlap,
    // nonlocal, rt-TDDFT hybrid gauge, local Pulay term and DeltaSpin).
    this->cal_operator_fs(ucell, gd, pv, pelec, dmat, psi, two_center_bundle,
                          orb, kv, isforce, isstress, cfg, td_stype, p_hamilt, parts, sparts);

    // MPI reduction for forces
    if (isforce)
    {
        Parallel_Reduce::reduce_pool(parts.fvl_dphi.c, parts.fvl_dphi.nr * parts.fvl_dphi.nc);
    }

    // MPI reduction for stresses
    if (isstress)
    {
        Parallel_Reduce::reduce_pool(sparts.svl_dphi.c, sparts.svl_dphi.nr * sparts.svl_dphi.nc);
    }

    // Handle DeePKS forces if enabled
    LCAO_domain::cal_deepks_fs(ucell, gd, pv, orb, kv, isforce, isstress, deepks, parts, sparts);

    // vdW force/stress and external-field forces
    LCAO_domain::cal_vdw_fields_fs(vdw_result, ucell, solvent, rhopw, locpp,
                                   isforce, isstress, parts, sparts);

    // DFT+U force/stress
    LCAO_domain::cal_dftu_fs(ucell, gd, pv, orb, kv, dmat, two_center_bundle, dftu, isforce, isstress, parts, sparts);


    // NOTE: finish_ftable is no longer needed as we don't use ForceStressArrays for overlap/kinetic
    // if (!PARAM.globalv.gamma_only_local)
    // {
    //     this->edm_cal.finish_ftable(fsr);
    // }

    // EXX force/stress
    LCAO_domain::cal_exx_fs(ucell, isforce, isstress, exx_info, exx_nao, parts, sparts);
    //--------------------------------
    // begin calculate and output force
    //--------------------------------
    if (isforce)
    {
        LCAO_domain::assemble_print_force(ucell, istestf, vdw_result, exx_info, symm, deepks.dpks_out_type, parts,
                                          force_invalid_threshold_ev, fcs);
    } // end of force calculation
    //---------------------------------
    // begin calculate and output stress
    //---------------------------------
    if (isstress)
    {
        LCAO_domain::assemble_print_stress(ucell, istests, vdw_result, exx_info, symm, deepks.dpks_out_type, sparts, scs);
    } // end of stress calculation

    ModuleBase::timer::end("Force_Stress_LCAO", "getForceStress");
    return;
}

// Operator-based force/stress terms: kinetic, overlap, nonlocal,
// rt-TDDFT hybrid gauge, local-potential Pulay term, and DeltaSpin.
template <typename T>
void Force_Stress_LCAO<T>::cal_operator_fs(UnitCell& ucell,
                                             const Grid_Driver& gd,
                                             Parallel_Orbitals& pv,
                                             const elecstate::ElecState* pelec,
                                             module_dm::Setup_DM<T>& dmat,
                                             const psi::Psi<T>* psi,
                                             const TwoCenterBundle& two_center_bundle,
                                             const LCAO_Orbitals& orb,
                                             const K_Vectors& kv,
                                             const bool isforce,
                                             const bool isstress,
                                             const FSCalcConfig& cfg,
                                             const int td_stype,
                                             hamilt::Hamilt<T>* p_hamilt,
                                             LCAOForceParts& parts,
                                             LCAOStressParts& sparts)
{

    // Calculate forces and stresses using new operator-based methods
    // Step 1: Calculate Energy Density Matrix (EDM) for overlap force
    // EDM = Σ_k w_k * ε_k * |ψ_k><ψ_k|
    module_dm::DensityMatrix<T, double> edm = edm_cal.cal_edm(pelec, *psi, *dmat.dm, kv, pv,
                                                           cfg.nspin, cfg.nbands, ucell, *this->RA);

    // Step 2: Handle different spin cases
    if (cfg.nspin == 1 || cfg.nspin == 2)
    {
        // For nspin=1 or nspin=2, use double precision
        // Switch to spin channel 1 for DMR access
        if (cfg.nspin == 2)
        {
            dmat.dm->switch_dmr(1);
            edm.switch_dmr(1);
        }

        const hamilt::HContainer<double>* dmR = dmat.dm->get_dmr_ptr(1);
        const hamilt::HContainer<double>* edmR = edm.get_dmr_ptr(1);

        // Calculate kinetic force/stress (uses DM)
        if (cfg.t_in_h)
        {
            hamilt::EKinetic<hamilt::OperatorLCAO<T, double>> tmp_ekinetic(
                nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
                two_center_bundle.kinetic_orb.get());
            tmp_ekinetic.cal_force_stress(isforce, isstress, dmR, parts.ftvnl_dphi, sparts.stvnl_dphi);
        }

        // Calculate overlap force/stress (uses EDM)
        hamilt::Overlap<hamilt::OperatorLCAO<T, double>> tmp_overlap(
            nullptr, kv.kvec_d, nullptr, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb.get());
        if(td_stype != 2)
        {
            tmp_overlap.cal_force_stress(isforce, isstress, edmR, parts.foverlap, sparts.soverlap);
        }

        // Calculate nonlocal force/stress (uses DM)
        hamilt::Nonlocal<hamilt::OperatorLCAO<T, double>> tmp_nonlocal(
            nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb_beta.get());
        tmp_nonlocal.cal_force_stress(isforce, isstress, dmR, parts.fvnl_dbeta, sparts.svnl_dbeta);
        
        if(td_stype == 2)
        {
            hamilt::TD_pot_hybrid<hamilt::OperatorLCAO<T, double>> tmp_hybrid(
                nullptr, &kv, nullptr, nullptr, orb, &ucell, orb.cutoffs(), &gd, nullptr);
            tmp_hybrid.cal_force_stress(isforce, dmR, parts.fpothybrid);

            cal_foverlap_rt(parts.foverlap, dmat, p_hamilt, kv, pv, ucell);
        }

        // Switch back to spin channel 0
        if (cfg.nspin == 2)
        {
            dmat.dm->switch_dmr(0);
            edm.switch_dmr(0);
        }

        // Calculate local potential force/stress (vl_dphi)
        // This uses grid integration, not operator-based method
        edm_cal.ParaV = &pv;
        PulayForceStress::cal_pulay_fs(parts.fvl_dphi, sparts.svl_dphi, *dmat.dm, ucell, pelec->pot,
                                       isforce, isstress, false /*reset dm to gint*/);
    }
    else if (cfg.nspin == 4)
    {

        // Kinetic force/stress from the complex DMR (nspin=4)
        if (cfg.t_in_h)
        {
            hamilt::EKinetic<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_ekinetic(
                nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
                two_center_bundle.kinetic_orb.get());
            tmp_ekinetic.cal_force_stress(isforce, isstress, dmat.dm->get_dmr_ptr(1), parts.ftvnl_dphi,
                                          sparts.stvnl_dphi);
        }

        // Overlap force/stress from the complex EDM (nspin=4)
        hamilt::Overlap<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_overlap(
            nullptr, kv.kvec_d, nullptr, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb.get());
        tmp_overlap.cal_force_stress(isforce, isstress, edm.get_dmr_ptr(1), parts.foverlap, sparts.soverlap);

        // For nspin=4 (non-collinear), need complex DMR
        // Create temporary complex DMR for DM
        hamilt::HContainer<std::complex<double>> tmp_dmr(dmat.dm->get_dmr_ptr(1)->get_paraV());
        std::vector<int> ijrs = dmat.dm->get_dmr_ptr(1)->get_ijr_info();
        tmp_dmr.insert_ijrs(&ijrs);
        tmp_dmr.allocate();
        dmat.dm->cal_dmr_full(&tmp_dmr, -1);
        // Nonlocal force/stress from the temporary complex DMR
        hamilt::Nonlocal<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_nonlocal(
            nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb_beta.get());
        tmp_nonlocal.cal_force_stress(isforce, isstress, &tmp_dmr, parts.fvnl_dbeta, sparts.svnl_dbeta);

        // Local-potential (vl_dphi) Pulay term via grid integration
        edm_cal.ParaV = &pv;
        PulayForceStress::cal_pulay_fs(parts.fvl_dphi, sparts.svl_dphi, *dmat.dm, ucell, pelec->pot,
                                       isforce, isstress, false);
    }

    // atomic force and stress for DeltaSpin
    if (cfg.sc_mag_switch)
    {
        if (isforce)
        {
            parts.force_dspin.create(ucell.nat, 3);
        }
        if (isstress)
        {
            sparts.stress_dspin.create(3, 3);
        }

        hamilt::DeltaSpin<hamilt::OperatorLCAO<T, double>> tmp_dspin(nullptr,
                                                                     kv.kvec_d,
                                                                     nullptr,
                                                                     ucell,
                                                                     &gd,
                                                                     two_center_bundle.overlap_orb_onsite.get(),
                                                                     orb.cutoffs());

        if (cfg.nspin == 2)
        {
            dmat.dm->switch_dmr(2);
        }
        const hamilt::HContainer<double>* dmr = dmat.dm->get_dmr_ptr(1);
        tmp_dspin.cal_force_stress(isforce, isstress, dmr, parts.force_dspin, sparts.stress_dspin);
        if (cfg.nspin == 2)
        {
            dmat.dm->switch_dmr(0);
        }
    }
}


#include "source_base/mathzone.h"

// local pseudopotential, ewald, core correction, scc terms in force
template <typename T>
void Force_Stress_LCAO<T>::calForcePwPart(UnitCell& ucell,
                                          ModuleBase::matrix& fvl_dvl,
                                          ModuleBase::matrix& fewalds,
                                          ModuleBase::matrix& fcc,
                                          ModuleBase::matrix& fscc,
                                          const double& etxc,
                                          const ModuleBase::matrix& vnew,
                                          const bool vnew_exist,
                                          const Charge* const chr,
                                          ModulePW::PW_Basis* rhopw,
                                          const pseudopot_cell_vl& locpp,
                                          const Structure_Factor& sf,
                                          const FSCalcConfig& cfg)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "calForcePwPart");
#ifdef __CUDA
    if (cfg.device == "gpu")
    {
        Forces<double, base_device::DEVICE_GPU> f_pw(nat);
        f_pw.cal_force_loc(ucell, fvl_dvl, rhopw, locpp.vloc, chr);
        f_pw.cal_force_ew(ucell, fewalds, rhopw, &sf);
        f_pw.cal_force_cc(fcc, rhopw, chr, locpp.numeric, ucell,
                          cfg.nspin, cfg.domag, cfg.domag_z, cfg.gga_grad);
        f_pw.cal_force_scc(fscc, rhopw, vnew, vnew_exist, locpp.numeric, ucell);
    }
    else
#endif
    {
        Forces<double, base_device::DEVICE_CPU> f_pw(nat);
        f_pw.cal_force_loc(ucell, fvl_dvl, rhopw, locpp.vloc, chr);
        f_pw.cal_force_ew(ucell, fewalds, rhopw, &sf);
        f_pw.cal_force_cc(fcc, rhopw, chr, locpp.numeric, ucell,
                          cfg.nspin, cfg.domag, cfg.domag_z, cfg.gga_grad);
        f_pw.cal_force_scc(fscc, rhopw, vnew, vnew_exist, locpp.numeric, ucell);
    }

    return;
}


template class Force_Stress_LCAO<double>;
template class Force_Stress_LCAO<std::complex<double>>;
