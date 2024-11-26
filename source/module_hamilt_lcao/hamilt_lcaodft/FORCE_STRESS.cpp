#include "FORCE_STRESS.h"

#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_dftu/dftu.h" //Quxin add for DFT+U on 20201029
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"
// new
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/potentials/efield.h"           // liuyu add 2022-05-18
#include "module_elecstate/potentials/gatefield.h"        // liuyu add 2022-09-13
#include "module_hamilt_general/module_surchem/surchem.h" //sunml add 2022-08-10
#include "module_hamilt_general/module_vdw/vdw.h"
#ifdef __DEEPKS
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h" //caoyu add for deepks 2021-06-03
#include "module_hamilt_lcao/module_deepks/LCAO_deepks_io.h" // mohan add 2024-07-22 
#endif
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dftu_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dspin_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/nonlocal_new.h"
#include "module_elecstate/elecstate_lcao.h"

template <typename T>
Force_Stress_LCAO<T>::Force_Stress_LCAO(Record_adj& ra, const int nat_in) : RA(&ra), f_pw(nat_in), nat(nat_in)
{
}
template <typename T>
Force_Stress_LCAO<T>::~Force_Stress_LCAO()
{
}
template <typename T>
void Force_Stress_LCAO<T>::getForceStress(const bool isforce,
                                          const bool isstress,
                                          const bool istestf,
                                          const bool istests,
                                          Parallel_Orbitals& pv,
                                          const elecstate::ElecState* pelec,
                                          const psi::Psi<T>* psi,
                                          Gint_Gamma& gint_gamma, // mohan add 2024-04-01
                                          Gint_k& gint_k,         // mohan add 2024-04-01
                                          const TwoCenterBundle& two_center_bundle,
                                          const LCAO_Orbitals& orb,
                                          ModuleBase::matrix& fcs,
                                          ModuleBase::matrix& scs,
                                          const Structure_Factor& sf,
                                          const K_Vectors& kv,
                                          ModulePW::PW_Basis* rhopw,
#ifdef __EXX
                                          Exx_LRI<double>& exx_lri_double,
                                          Exx_LRI<std::complex<double>>& exx_lri_complex,
#endif
                                          ModuleSymmetry::Symmetry* symm)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "getForceStress");
    ModuleBase::timer::tick("Force_Stress_LCAO", "getForceStress");

    if (!isforce && !isstress)
    {
        ModuleBase::timer::tick("Force_Stress_LCAO", "getForceStress");
        return;
    }

    const int nat = GlobalC::ucell.nat;

    ForceStressArrays fsr; // mohan add 2024-06-15

    // total force : ModuleBase::matrix fcs;

    // part of total force
    ModuleBase::matrix foverlap;
    ModuleBase::matrix ftvnl_dphi;
    ModuleBase::matrix fvnl_dbeta;
    ModuleBase::matrix fvl_dphi;
    ModuleBase::matrix fvl_dvl;
    ModuleBase::matrix fewalds;
    ModuleBase::matrix fcc;
    ModuleBase::matrix fscc;

    fvl_dphi.create(nat, 3); // must do it now, update it later, noted by zhengdy

    if (isforce)
    {
        fcs.create(nat, 3);
        foverlap.create(nat, 3);
        ftvnl_dphi.create(nat, 3);
        fvnl_dbeta.create(nat, 3);
        fvl_dvl.create(nat, 3);
        fewalds.create(nat, 3);
        fcc.create(nat, 3);
        fscc.create(nat, 3);

        // calculate basic terms in Force, same method with PW base
        this->calForcePwPart(fvl_dvl,
                             fewalds,
                             fcc,
                             fscc,
                             pelec->f_en.etxc,
                             pelec->vnew,
                             pelec->vnew_exist,
                             pelec->charge,
                             rhopw,
                             sf);
    }

    // total stress : ModuleBase::matrix scs
    ModuleBase::matrix sigmacc;
    ModuleBase::matrix sigmadvl;
    ModuleBase::matrix sigmaewa;
    ModuleBase::matrix sigmaxc;
    ModuleBase::matrix sigmahar;
    ModuleBase::matrix soverlap;
    ModuleBase::matrix stvnl_dphi;
    ModuleBase::matrix svnl_dbeta;
    ModuleBase::matrix svl_dphi;
#ifdef __DEEPKS
    ModuleBase::matrix svnl_dalpha; // deepks
#endif

    //! stress
    if (isstress)
    {
        scs.create(3, 3);
        sigmacc.create(3, 3);
        sigmadvl.create(3, 3);
        sigmaewa.create(3, 3);
        sigmaxc.create(3, 3);
        sigmahar.create(3, 3);

        soverlap.create(3, 3);
        stvnl_dphi.create(3, 3);
        svnl_dbeta.create(3, 3);
        svl_dphi.create(3, 3);
#ifdef __DEEPKS
        svnl_dalpha.create(3, 3);
#endif
        // calculate basic terms in Stress, similar method with PW base
        this->calStressPwPart(sigmadvl,
                              sigmahar,
                              sigmaewa,
                              sigmacc,
                              sigmaxc,
                              pelec->f_en.etxc,
                              pelec->charge,
                              rhopw,
                              sf);
    }

    //! atomic forces from integration (4 terms)
    this->integral_part(PARAM.globalv.gamma_only_local,
                        isforce,
                        isstress,
                        fsr,
                        pelec,
                        psi,
                        foverlap,
                        ftvnl_dphi,
                        fvnl_dbeta,
                        fvl_dphi,
                        soverlap,
                        stvnl_dphi,
                        svnl_dbeta,
                        svl_dphi,
#ifdef __DEEPKS
                        svnl_dalpha,
#endif
                        gint_gamma,
                        gint_k,
                        two_center_bundle,
                        orb,
                        pv,
                        kv);
    // calculate force and stress for Nonlocal part
    if(PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
    {
        hamilt::NonlocalNew<hamilt::OperatorLCAO<T, double>> tmp_nonlocal(
                    nullptr,
                    kv.kvec_d,
                    nullptr,
                    &GlobalC::ucell,
                    orb.cutoffs(),
                    &GlobalC::GridD,
                    two_center_bundle.overlap_orb_beta.get()
            );

        const auto* dm_p = dynamic_cast<const elecstate::ElecStateLCAO<T>*>(pelec)->get_DM();
        if(PARAM.inp.nspin == 2)
        {
            const_cast<elecstate::DensityMatrix<T, double>*>(dm_p)->switch_dmr(1);
        }
        const hamilt::HContainer<double>* dmr = dm_p->get_DMR_pointer(1);
        tmp_nonlocal.cal_force_stress(isforce, isstress, dmr, fvnl_dbeta, svnl_dbeta);
        if(PARAM.inp.nspin == 2)
        {
            const_cast<elecstate::DensityMatrix<T, double>*>(dm_p)->switch_dmr(0);
        }
    }
    else if(PARAM.inp.nspin == 4)
    {
        hamilt::NonlocalNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_nonlocal(
                    nullptr,
                    kv.kvec_d,
                    nullptr,
                    &GlobalC::ucell,
                    orb.cutoffs(),
                    &GlobalC::GridD,
                    two_center_bundle.overlap_orb_beta.get()
            );

        // calculate temporary complex DMR for nonlocal force&stress
        // In fact, only SOC part need the imaginary part of DMR for correct force&stress
        const auto* dm_p = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)->get_DM();
        hamilt::HContainer<std::complex<double>> tmp_dmr(dm_p->get_DMR_pointer(1)->get_paraV());
        std::vector<int> ijrs = dm_p->get_DMR_pointer(1)->get_ijr_info();
        tmp_dmr.insert_ijrs(&ijrs);
        tmp_dmr.allocate();
        dm_p->cal_DMR_full(&tmp_dmr);
        tmp_nonlocal.cal_force_stress(isforce, isstress, &tmp_dmr, fvnl_dbeta, svnl_dbeta);
    }
    

    //! forces and stress from vdw
    //  Peize Lin add 2014-04-04, update 2021-03-09
    //  jiyy add 2019-05-18, update 2021-05-02
    ModuleBase::matrix force_vdw;
    ModuleBase::matrix stress_vdw;
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, PARAM.inp);
    if (vdw_solver != nullptr)
    {
        if (isforce)
        {
            force_vdw.create(nat, 3);
            const std::vector<ModuleBase::Vector3<double>>& force_vdw_temp = vdw_solver->get_force();
            for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
            {
                force_vdw(iat, 0) = force_vdw_temp[iat].x;
                force_vdw(iat, 1) = force_vdw_temp[iat].y;
                force_vdw(iat, 2) = force_vdw_temp[iat].z;
            }
        }
        if (isstress)
        {
            stress_vdw = vdw_solver->get_stress().to_matrix();
        }
    }

    //! forces from E-field
    ModuleBase::matrix fefield;
    if (PARAM.inp.efield_flag && isforce)
    {
        fefield.create(nat, 3);
        elecstate::Efield::compute_force(GlobalC::ucell, fefield);
    }

    //! atomic forces from E-field of rt-TDDFT
    ModuleBase::matrix fefield_tddft;
    if (PARAM.inp.esolver_type == "TDDFT" && isforce)
    {
        fefield_tddft.create(nat, 3);
        elecstate::Efield::compute_force(GlobalC::ucell, fefield_tddft);
    }

    //! atomic forces from gate field
    ModuleBase::matrix fgate;
    if (PARAM.inp.gate_flag && isforce)
    {
        fgate.create(nat, 3);
        elecstate::Gatefield::compute_force(GlobalC::ucell, fgate);
    }

    //! atomic forces from implicit solvation model
    ModuleBase::matrix fsol;
    if (PARAM.inp.imp_sol && isforce)
    {
        fsol.create(nat, 3);
        GlobalC::solvent_model.cal_force_sol(GlobalC::ucell, rhopw, fsol);
    }

    //! atomic forces from DFT+U (Quxin version)
    ModuleBase::matrix force_dftu;
    ModuleBase::matrix stress_dftu;

    if (PARAM.inp.dft_plus_u) // Quxin add for DFT+U on 20201029
    {
        if (isforce)
        {
            force_dftu.create(nat, 3);
        }
        if (isstress)
        {
            stress_dftu.create(3, 3);
        }
        if (PARAM.inp.dft_plus_u == 2)
        {
            GlobalC::dftu.force_stress(pelec,
                                       pv,
                                       fsr, // mohan 2024-06-16
                                       force_dftu,
                                       stress_dftu,
                                       kv);
        }
        else
        {
            hamilt::DFTU<hamilt::OperatorLCAO<T, double>> tmp_dftu(nullptr, // HK and SK are not used for force&stress
                                                                   kv.kvec_d,
                                                                   nullptr, // HR are not used for force&stress
                                                                   GlobalC::ucell,
                                                                   &GlobalC::GridD,
                                                                   two_center_bundle.overlap_orb_onsite.get(),
                                                                   orb.cutoffs(),
                                                                   &GlobalC::dftu);

            tmp_dftu.cal_force_stress(isforce, isstress, force_dftu, stress_dftu);
        }
    }

    // atomic force and stress for DeltaSpin
    ModuleBase::matrix force_dspin;
    ModuleBase::matrix stress_dspin;
    if(PARAM.inp.sc_mag_switch)
    {
        if (isforce)
        {
            force_dspin.create(nat, 3);
        }
        if (isstress)
        {
            stress_dspin.create(3, 3);
        }

        hamilt::DeltaSpin<hamilt::OperatorLCAO<T, double>> tmp_dspin(
                    nullptr,
                    kv.kvec_d,
                    nullptr,
                    GlobalC::ucell,
                    &GlobalC::GridD,
                    two_center_bundle.overlap_orb_onsite.get(),
                    orb.cutoffs()
            );

        const auto* dm_p = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)->get_DM();
        if(PARAM.inp.nspin == 2)
        {
            const_cast<elecstate::DensityMatrix<std::complex<double>, double>*>(dm_p)->switch_dmr(2);
        }
        const hamilt::HContainer<double>* dmr = dm_p->get_DMR_pointer(1);
        tmp_dspin.cal_force_stress(isforce, isstress, dmr, force_dspin, stress_dspin);
        if(PARAM.inp.nspin == 2)
        {
            const_cast<elecstate::DensityMatrix<std::complex<double>, double>*>(dm_p)->switch_dmr(0);
        }
    }

    if (!PARAM.globalv.gamma_only_local)
    {
        this->flk.finish_ftable(fsr);
    }

#ifdef __EXX
    // Force and Stress contribution from exx
    ModuleBase::matrix force_exx;
    ModuleBase::matrix stress_exx;
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (isforce)
        {
            if (GlobalC::exx_info.info_ri.real_number)
            {
                exx_lri_double.cal_exx_force();
                force_exx = GlobalC::exx_info.info_global.hybrid_alpha * exx_lri_double.force_exx;
            }
            else
            {
                exx_lri_complex.cal_exx_force();
                force_exx = GlobalC::exx_info.info_global.hybrid_alpha * exx_lri_complex.force_exx;
            }
        }
        if (isstress)
        {
            if (GlobalC::exx_info.info_ri.real_number)
            {
                exx_lri_double.cal_exx_stress();
                stress_exx = GlobalC::exx_info.info_global.hybrid_alpha * exx_lri_double.stress_exx;
            }
            else
            {
                exx_lri_complex.cal_exx_stress();
                stress_exx = GlobalC::exx_info.info_global.hybrid_alpha * exx_lri_complex.stress_exx;
            }
        }
    }
#endif
    //--------------------------------
    // begin calculate and output force
    //--------------------------------
    if (isforce)
    {
        //---------------------------------
        // sum all parts of force!
        //---------------------------------
        for (int i = 0; i < 3; i++)
        {
            double sum = 0.0;

            for (int iat = 0; iat < nat; iat++)
            {
                fcs(iat, i) += foverlap(iat, i) + ftvnl_dphi(iat, i) + fvnl_dbeta(iat, i) + fvl_dphi(iat, i)
                               + fvl_dvl(iat, i) // derivative of local potential force (pw)
                               + fewalds(iat, i) // ewald force (pw)
                               + fcc(iat, i)     // nonlinear core correction force (pw)
                               + fscc(iat, i);   // self consistent corretion force (pw)

                // Force contribution from DFT+U, Quxin add on 20201029
                if (PARAM.inp.dft_plus_u)
                {
                    fcs(iat, i) += force_dftu(iat, i);
                }
                if (PARAM.inp.sc_mag_switch)
                {
                    fcs(iat, i) += force_dspin(iat, i);
                }
#ifdef __EXX
                // Force contribution from exx
                if (GlobalC::exx_info.info_global.cal_exx)
                {
                    fcs(iat, i) += force_exx(iat, i);
                }
#endif
                // VDW force of vdwd2 or vdwd3
                if (vdw_solver != nullptr)
                {
                    fcs(iat, i) += force_vdw(iat, i);
                }
                // E-field force
                if (PARAM.inp.efield_flag)
                {
                    fcs(iat, i) += fefield(iat, i);
                }
                // E-field force of tddft
                if (PARAM.inp.esolver_type == "TDDFT")
                {
                    fcs(iat, i) += fefield_tddft(iat, i);
                }
                // Gate field force
                if (PARAM.inp.gate_flag)
                {
                    fcs(iat, i) += fgate(iat, i);
                }
                // implicit solvation model
                if (PARAM.inp.imp_sol)
                {
                    fcs(iat, i) += fsol(iat, i);
                }
#ifdef __DEEPKS
                // mohan add 2021-08-04
                if (PARAM.inp.deepks_scf)
                {
                    fcs(iat, i) += GlobalC::ld.F_delta(iat, i);
                }
#endif
                // sum total force for correction
                sum += fcs(iat, i);
            }

            if (!(PARAM.inp.gate_flag || PARAM.inp.efield_flag))
            {
                for (int iat = 0; iat < nat; ++iat)
                {
                    fcs(iat, i) -= sum / nat;
                }
            }

            // xiaohui add "OUT_LEVEL", 2015-09-16
            if (PARAM.inp.out_level != "m")
            {
                GlobalV::ofs_running << " correction force for each atom along direction " << i + 1 << " is "
                                     << sum / nat << std::endl;
            }
        }

        if (PARAM.inp.gate_flag || PARAM.inp.efield_flag)
        {
            GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
        }

        // pengfei 2016-12-20
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            this->forceSymmetry(fcs, symm);
        }

#ifdef __DEEPKS
        // DeePKS force
        if (PARAM.inp.deepks_out_labels) // not parallelized yet
        {
            const std::string file_ftot = PARAM.globalv.global_out_dir + "deepks_ftot.npy";
            LCAO_deepks_io::save_npy_f(fcs, 
                                       file_ftot, 
                                       GlobalC::ucell.nat, 
                                       GlobalV::MY_RANK); // Ty/Bohr, F_tot

            if (PARAM.inp.deepks_scf)
            {
                const std::string file_fbase = PARAM.globalv.global_out_dir + "deepks_fbase.npy";
                LCAO_deepks_io::save_npy_f(fcs - GlobalC::ld.F_delta, 
                                           file_fbase, 
                                           GlobalC::ucell.nat, 
                                           GlobalV::MY_RANK); // Ry/Bohr, F_base

                if (!PARAM.inp.deepks_equiv) // training with force label not supported by equivariant version now
                {
                    if (PARAM.globalv.gamma_only_local)
                    {
                        const std::vector<std::vector<double>>& dm_gamma
                            = dynamic_cast<const elecstate::ElecStateLCAO<double>*>(pelec)->get_DM()->get_DMK_vector();
                        GlobalC::ld.cal_gdmx(dm_gamma[0], GlobalC::ucell, orb, GlobalC::GridD, isstress);
                    }
                    else
                    {
                        const std::vector<std::vector<std::complex<double>>>& dm_k
                            = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)
                                  ->get_DM()
                                  ->get_DMK_vector();

                        GlobalC::ld.cal_gdmx_k(dm_k,
                                               GlobalC::ucell,
                                               orb,
                                               GlobalC::GridD,
                                               kv.get_nks(),
                                               kv.kvec_d,
                                               isstress);
                    }
                    if (PARAM.inp.deepks_out_unittest)
                    {
                        GlobalC::ld.check_gdmx(GlobalC::ucell.nat);
                    }
                    GlobalC::ld.cal_gvx(GlobalC::ucell.nat);

                    if (PARAM.inp.deepks_out_unittest)
                    {
                        GlobalC::ld.check_gvx(GlobalC::ucell.nat);
                    }

                    LCAO_deepks_io::save_npy_gvx(GlobalC::ucell.nat,
                      GlobalC::ld.des_per_atom,
                      GlobalC::ld.gvx_tensor,
                      PARAM.globalv.global_out_dir,
                      GlobalV::MY_RANK);
                }
            }
            else
            {
                const std::string file_fbase = PARAM.globalv.global_out_dir + "deepks_fbase.npy";
                LCAO_deepks_io::save_npy_f(fcs, 
                  file_fbase, 
                  GlobalC::ucell.nat,
                  GlobalV::MY_RANK); // no scf, F_base=F_tot
            }
        }
#endif
        // print Rydberg force or not
        bool ry = false;
        if (istestf)
        {
            // test
            // ModuleBase::matrix fvlocal;
            // fvlocal.create(nat,3);
            ModuleBase::matrix ftvnl;
            ftvnl.create(nat, 3);
            for (int iat = 0; iat < nat; iat++)
            {
                for (int i = 0; i < 3; i++)
                {
                    // fvlocal(iat,i) = fvl_dphi(iat,i) + fvl_dvl(iat,i);
                    ftvnl(iat, i) = ftvnl_dphi(iat, i) + fvnl_dbeta(iat, i);
                }
            }

            GlobalV::ofs_running << "\n PARTS OF FORCE: " << std::endl;
            GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
            GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(8) << std::endl;
            //-----------------------------
            // regular force terms test.
            //-----------------------------
            // this->print_force("OVERLAP    FORCE",foverlap,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "OVERLAP    FORCE", foverlap, false);
            //  this->print_force("TVNL_DPHI  force",ftvnl_dphi,PARAM.inp.test_force);
            //  this->print_force("VNL_DBETA  force",fvnl_dbeta,PARAM.inp.test_force);
            // this->print_force("T_VNL      FORCE",ftvnl,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "T_VNL      FORCE", ftvnl, false);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VL_dPHI    FORCE", fvl_dphi, false);
            // this->print_force("VL_dPHI    FORCE",fvl_dphi,1,ry);
            // this->print_force("VL_dVL     FORCE",fvl_dvl,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VL_dVL     FORCE", fvl_dvl, false);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EWALD      FORCE", fewalds, false);
            // 	this->print_force("VLOCAL     FORCE",fvlocal,PARAM.inp.test_force);
            // this->print_force("EWALD      FORCE",fewalds,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "NLCC       FORCE", fcc, false);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "SCC        FORCE", fscc, false);
            // this->print_force("NLCC       FORCE",fcc,1,ry);
            // this->print_force("SCC        FORCE",fscc,1,ry);
            //-------------------------------
            // put extra force here for test!
            //-------------------------------
            if (PARAM.inp.efield_flag)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EFIELD     FORCE", fefield, false);
                // this->print_force("EFIELD     FORCE",fefield,1,ry);
            }
            if (PARAM.inp.esolver_type == "TDDFT")
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EFIELD_TDDFT     FORCE", fefield_tddft, false);
                // this->print_force("EFIELD_TDDFT     FORCE",fefield_tddft,1,ry);
            }
            if (PARAM.inp.gate_flag)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "GATEFIELD     FORCE", fgate, false);
                // this->print_force("GATEFIELD     FORCE",fgate,1,ry);
            }
            if (PARAM.inp.imp_sol)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "IMP_SOL     FORCE", fsol, false);
                // this->print_force("IMP_SOL     FORCE",fsol,1,ry);
            }
            if (vdw_solver != nullptr)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VDW        FORCE", force_vdw, false);
                // this->print_force("VDW        FORCE",force_vdw,1,ry);
            }
            if (PARAM.inp.dft_plus_u)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "DFT+U      FORCE", force_dftu, false);
            }
            if (PARAM.inp.sc_mag_switch)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "DeltaSpin  FORCE", force_dspin, false);
            }
#ifdef __DEEPKS
            // caoyu add 2021-06-03
            if (PARAM.inp.deepks_scf)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "DeePKS 	FORCE", GlobalC::ld.F_delta, true);
                // this->print_force("DeePKS 	FORCE", GlobalC::ld.F_delta, 1, ry);
            }
#endif
        }

        GlobalV::ofs_running << std::setiosflags(std::ios::left);

        // this->printforce_total(ry, istestf, fcs);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "TOTAL-FORCE (eV/Angstrom)", fcs, false);
        if (istestf)
        {
            GlobalV::ofs_running << "\n FORCE INVALID TABLE." << std::endl;
            GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(5) << "x" << std::setw(5) << "y"
                                 << std::setw(5) << "z" << std::endl;
            for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
            {
                GlobalV::ofs_running << " " << std::setw(8) << iat;
                for (int i = 0; i < 3; i++)
                {
                    if (std::abs(fcs(iat, i) * ModuleBase::Ry_to_eV / 0.529177)
                        < Force_Stress_LCAO::force_invalid_threshold_ev)
                    {
                        fcs(iat, i) = 0.0;
                        GlobalV::ofs_running << std::setw(5) << "1";
                    }
                    else
                    {
                        GlobalV::ofs_running << std::setw(5) << "0";
                    }
                }
                GlobalV::ofs_running << std::endl;
            }
        }
    } // end of force calculation
    //---------------------------------
    // begin calculate and output stress
    //---------------------------------
    if (isstress)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                scs(i, j) += soverlap(i, j) + stvnl_dphi(i, j) + svnl_dbeta(i, j) + svl_dphi(i, j)
                             + sigmadvl(i, j)  // derivative of local potential stress (pw)
                             + sigmaewa(i, j)  // ewald stress (pw)
                             + sigmacc(i, j)   // nonlinear core correction stress (pw)
                             + sigmaxc(i, j)   // exchange corretion stress
                             + sigmahar(i, j); // hartree stress

                // VDW stress from linpz and jiyy
                if (vdw_solver != nullptr)
                {
                    scs(i, j) += stress_vdw(i, j);
                }
                // DFT plus U stress from qux
                if (PARAM.inp.dft_plus_u)
                {
                    scs(i, j) += stress_dftu(i, j);
                }
                if (PARAM.inp.sc_mag_switch)
                {
                    scs(i, j) += stress_dspin(i, j);
                }
#ifdef __EXX
                // Stress contribution from exx
                if (GlobalC::exx_info.info_global.cal_exx)
                {
                    scs(i, j) += stress_exx(i, j);
                }
#endif
            }
        }
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            symm->symmetrize_mat3(scs, GlobalC::ucell.lat);
        } // end symmetry

#ifdef __DEEPKS
        if (PARAM.inp.deepks_out_labels) // not parallelized yet
        {
            const std::string file_s = PARAM.globalv.global_out_dir + "deepks_sbase.npy";
            LCAO_deepks_io::save_npy_s(scs,
                                      file_s,
                                      GlobalC::ucell.omega,
                                      GlobalV::MY_RANK); // change to energy unit Ry when printing, S_base;
        }
        if (PARAM.inp.deepks_scf)
        {
            if (ModuleSymmetry::Symmetry::symm_flag == 1)
            {
                symm->symmetrize_mat3(svnl_dalpha, GlobalC::ucell.lat);
            } // end symmetry
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    scs(i, j) += svnl_dalpha(i, j);
                }
            }
        }
        if (PARAM.inp.deepks_out_labels) // not parallelized yet
        {
			const std::string file_s = PARAM.globalv.global_out_dir + "deepks_stot.npy";
			LCAO_deepks_io::save_npy_s(
					scs,
					file_s,
					GlobalC::ucell.omega,
					GlobalV::MY_RANK); // change to energy unit Ry when printing, S_tot, w/ model

            // wenfei add 2021/11/2
            if (PARAM.inp.deepks_scf)
            {

                if (!PARAM.inp.deepks_equiv) // training with stress label not supported by equivariant version now
                {
                    GlobalC::ld.cal_gvepsl(GlobalC::ucell.nat);

                    LCAO_deepks_io::save_npy_gvepsl(
                      GlobalC::ucell.nat, 
                      GlobalC::ld.des_per_atom,
                      GlobalC::ld.gvepsl_tensor,
                      PARAM.globalv.global_out_dir,
                      GlobalV::MY_RANK); //  unitless, grad_vepsl
                }
            }
        }
#endif

        // print Rydberg stress or not
        bool ry = false;

        // test stress each terms if needed
        if (istests)
        {
            // test
            ModuleBase::matrix svlocal;
            svlocal.create(3, 3);
            ModuleBase::matrix stvnl;
            stvnl.create(3, 3);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    svlocal(i, j) = svl_dphi(i, j) + sigmadvl(i, j);
                    stvnl(i, j) = stvnl_dphi(i, j) + svnl_dbeta(i, j);
                }
            }

            GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
            GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
            GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(10) << std::endl;
            ModuleIO::print_stress("OVERLAP  STRESS", soverlap, PARAM.inp.test_stress, ry);
            // test
            ModuleIO::print_stress("T        STRESS", stvnl_dphi, PARAM.inp.test_stress, ry);
            ModuleIO::print_stress("VNL      STRESS", svnl_dbeta, PARAM.inp.test_stress, ry);

            ModuleIO::print_stress("T_VNL    STRESS", stvnl, PARAM.inp.test_stress, ry);

            ModuleIO::print_stress("VL_dPHI  STRESS", svl_dphi, PARAM.inp.test_stress, ry);
            ModuleIO::print_stress("VL_dVL   STRESS", sigmadvl, PARAM.inp.test_stress, ry);
            ModuleIO::print_stress("HAR      STRESS", sigmahar, PARAM.inp.test_stress, ry);

            ModuleIO::print_stress("EWALD    STRESS", sigmaewa, PARAM.inp.test_stress, ry);
            ModuleIO::print_stress("cc       STRESS", sigmacc, PARAM.inp.test_stress, ry);
            //		ModuleIO::print_stress("NLCC       STRESS",sigmacc,PARAM.inp.test_stress,ry);
            ModuleIO::print_stress("XC       STRESS", sigmaxc, PARAM.inp.test_stress, ry);
            if (vdw_solver != nullptr)
            {
                ModuleIO::print_stress("VDW      STRESS", sigmaxc, PARAM.inp.test_stress, ry);
            }
            if (PARAM.inp.dft_plus_u)
            {
                ModuleIO::print_stress("DFTU     STRESS", stress_dftu, PARAM.inp.test_stress, ry);
            }
            if (PARAM.inp.sc_mag_switch)
            {
                ModuleIO::print_stress("DeltaSpin  STRESS", stress_dspin, PARAM.inp.test_stress, ry);
            }
            ModuleIO::print_stress("TOTAL    STRESS", scs, PARAM.inp.test_stress, ry);

        } // end of test
        GlobalV::ofs_running << std::setiosflags(std::ios::left);
        // print total stress
        ModuleIO::print_stress("TOTAL-STRESS", scs, true, ry);

        double unit_transform = 0.0;
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
        double external_stress[3] = {PARAM.inp.press1, PARAM.inp.press2, PARAM.inp.press3};

        for (int i = 0; i < 3; i++)
        {
            scs(i, i) -= external_stress[i] / unit_transform;
        }
    } // end of stress calculation

    ModuleBase::timer::tick("Force_Stress_LCAO", "getForceStress");
    return;
}

// local pseudopotential, ewald, core correction, scc terms in force
template <typename T>
void Force_Stress_LCAO<T>::calForcePwPart(ModuleBase::matrix& fvl_dvl,
                                          ModuleBase::matrix& fewalds,
                                          ModuleBase::matrix& fcc,
                                          ModuleBase::matrix& fscc,
                                          const double& etxc,
                                          const ModuleBase::matrix& vnew,
                                          const bool vnew_exist,
                                          const Charge* const chr,
                                          ModulePW::PW_Basis* rhopw,
                                          const Structure_Factor& sf)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "calForcePwPart");
    //--------------------------------------------------------
    // local pseudopotential force:
    // use charge density; plane wave; local pseudopotential;
    //--------------------------------------------------------
    f_pw.cal_force_loc(fvl_dvl, rhopw, chr);
    //--------------------------------------------------------
    // ewald force: use plane wave only.
    //--------------------------------------------------------
    f_pw.cal_force_ew(fewalds, rhopw, &sf); // remain problem

    //--------------------------------------------------------
    // force due to core correlation.
    //--------------------------------------------------------
    f_pw.cal_force_cc(fcc, rhopw, chr,GlobalC::ucell);
    //--------------------------------------------------------
    // force due to self-consistent charge.
    //--------------------------------------------------------
    f_pw.cal_force_scc(fscc, rhopw, vnew, vnew_exist,GlobalC::ucell);
    return;
}

// overlap, kinetic, nonlocal pseudopotential, Local potential terms in force and stress
template <>
void Force_Stress_LCAO<double>::integral_part(const bool isGammaOnly,
                                              const bool isforce,
                                              const bool isstress,
                                              ForceStressArrays& fsr, // mohan add 2024-06-15
                                              const elecstate::ElecState* pelec,
                                              const psi::Psi<double>* psi,
                                              ModuleBase::matrix& foverlap,
                                              ModuleBase::matrix& ftvnl_dphi,
                                              ModuleBase::matrix& fvnl_dbeta,
                                              ModuleBase::matrix& fvl_dphi,
                                              ModuleBase::matrix& soverlap,
                                              ModuleBase::matrix& stvnl_dphi,
                                              ModuleBase::matrix& svnl_dbeta,
                                              ModuleBase::matrix& svl_dphi,
#if __DEEPKS
                                              ModuleBase::matrix& svnl_dalpha,
#endif
                                              Gint_Gamma& gint_gamma, // mohan add 2024-04-01
                                              Gint_k& gint_k,         // mohan add 2024-04-01
                                              const TwoCenterBundle& two_center_bundle,
                                              const LCAO_Orbitals& orb,
                                              const Parallel_Orbitals& pv,
                                              const K_Vectors& kv)
{

    flk.ftable(isforce,
               isstress,
               fsr, // mohan add 2024-06-15
               GlobalC::ucell,
               psi,
               pelec,
               foverlap,
               ftvnl_dphi,
               fvnl_dbeta,
               fvl_dphi,
               soverlap,
               stvnl_dphi,
               svnl_dbeta,
               svl_dphi,
#if __DEEPKS
               svnl_dalpha,
#endif
               gint_gamma,
               two_center_bundle,
               orb,
               pv);
    return;
}

template <>
void Force_Stress_LCAO<std::complex<double>>::integral_part(const bool isGammaOnly,
                                                            const bool isforce,
                                                            const bool isstress,
                                                            ForceStressArrays& fsr, // mohan add 2024-06-15
                                                            const elecstate::ElecState* pelec,
                                                            const psi::Psi<std::complex<double>>* psi,
                                                            ModuleBase::matrix& foverlap,
                                                            ModuleBase::matrix& ftvnl_dphi,
                                                            ModuleBase::matrix& fvnl_dbeta,
                                                            ModuleBase::matrix& fvl_dphi,
                                                            ModuleBase::matrix& soverlap,
                                                            ModuleBase::matrix& stvnl_dphi,
                                                            ModuleBase::matrix& svnl_dbeta,
                                                            ModuleBase::matrix& svl_dphi,
#if __DEEPKS
                                                            ModuleBase::matrix& svnl_dalpha,
#endif
                                                            Gint_Gamma& gint_gamma,
                                                            Gint_k& gint_k,
                                                            const TwoCenterBundle& two_center_bundle,
                                                            const LCAO_Orbitals& orb,
                                                            const Parallel_Orbitals& pv,
                                                            const K_Vectors& kv)
{
    flk.ftable(isforce,
               isstress,
               fsr, // mohan add 2024-06-16
               GlobalC::ucell,
               psi,
               pelec,
               foverlap,
               ftvnl_dphi,
               fvnl_dbeta,
               fvl_dphi,
               soverlap,
               stvnl_dphi,
               svnl_dbeta,
               svl_dphi,
#if __DEEPKS
               svnl_dalpha,
#endif
               gint_k,
               two_center_bundle,
               orb,
               pv,
               &kv,
               this->RA);
    return;
}

// vlocal, hartree, ewald, core correction, exchange-correlation terms in stress
template <typename T>
void Force_Stress_LCAO<T>::calStressPwPart(ModuleBase::matrix& sigmadvl,
                                           ModuleBase::matrix& sigmahar,
                                           ModuleBase::matrix& sigmaewa,
                                           ModuleBase::matrix& sigmacc,
                                           ModuleBase::matrix& sigmaxc,
                                           const double& etxc,
                                           const Charge* const chr,
                                           ModulePW::PW_Basis* rhopw,
                                           const Structure_Factor& sf)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "calStressPwPart");
    //--------------------------------------------------------
    // local pseudopotential stress:
    // use charge density; plane wave; local pseudopotential;
    //--------------------------------------------------------
    sc_pw.stress_loc(sigmadvl, rhopw, &sf, 0, chr);

    //--------------------------------------------------------
    // hartree term
    //--------------------------------------------------------
    sc_pw.stress_har(sigmahar, rhopw, 0, chr);

    //--------------------------------------------------------
    // ewald stress: use plane wave only.
    //--------------------------------------------------------
    sc_pw.stress_ewa(sigmaewa, rhopw, 0); // remain problem

    //--------------------------------------------------------
    // stress due to core correlation.
    //--------------------------------------------------------
    sc_pw.stress_cc(sigmacc, rhopw, &sf, 0, chr);

    //--------------------------------------------------------
    // stress due to self-consistent charge.
    //--------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
        sigmaxc(i, i) = -etxc / GlobalC::ucell.omega;
    }
    // Exchange-correlation for PBE
    sc_pw.stress_gga(sigmaxc, rhopw, chr);

    return;
}

#include "module_base/mathzone.h"
// do symmetry for total force
template <typename T>
void Force_Stress_LCAO<T>::forceSymmetry(ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm)
{
    double d1, d2, d3;
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        ModuleBase::Mathzone::Cartesian_to_Direct(fcs(iat, 0),
                                                  fcs(iat, 1),
                                                  fcs(iat, 2),
                                                  GlobalC::ucell.a1.x,
                                                  GlobalC::ucell.a1.y,
                                                  GlobalC::ucell.a1.z,
                                                  GlobalC::ucell.a2.x,
                                                  GlobalC::ucell.a2.y,
                                                  GlobalC::ucell.a2.z,
                                                  GlobalC::ucell.a3.x,
                                                  GlobalC::ucell.a3.y,
                                                  GlobalC::ucell.a3.z,
                                                  d1,
                                                  d2,
                                                  d3);

        fcs(iat, 0) = d1;
        fcs(iat, 1) = d2;
        fcs(iat, 2) = d3;
    }
    symm->symmetrize_vec3_nat(fcs.c);
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        ModuleBase::Mathzone::Direct_to_Cartesian(fcs(iat, 0),
                                                  fcs(iat, 1),
                                                  fcs(iat, 2),
                                                  GlobalC::ucell.a1.x,
                                                  GlobalC::ucell.a1.y,
                                                  GlobalC::ucell.a1.z,
                                                  GlobalC::ucell.a2.x,
                                                  GlobalC::ucell.a2.y,
                                                  GlobalC::ucell.a2.z,
                                                  GlobalC::ucell.a3.x,
                                                  GlobalC::ucell.a3.y,
                                                  GlobalC::ucell.a3.z,
                                                  d1,
                                                  d2,
                                                  d3);

        fcs(iat, 0) = d1;
        fcs(iat, 1) = d2;
        fcs(iat, 2) = d3;
    }
    return;
}

template class Force_Stress_LCAO<double>;
template class Force_Stress_LCAO<std::complex<double>>;
