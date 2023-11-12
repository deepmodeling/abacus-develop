#include "FORCE_STRESS.h"

#include "module_hamilt_lcao/module_dftu/dftu.h" //Quxin add for DFT+U on 20201029
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"
// new
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/potentials/efield.h"           // liuyu add 2022-05-18
#include "module_elecstate/potentials/gatefield.h"        // liuyu add 2022-09-13
#include "module_hamilt_general/module_surchem/surchem.h" //sunml add 2022-08-10
#include "module_hamilt_general/module_vdw/vdw.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h" //caoyu add for deepks 2021-06-03
#include "module_elecstate/elecstate_lcao.h"
#endif

template<typename T>
double Force_Stress_LCAO<T>::force_invalid_threshold_ev = 0.00;
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
                                          Local_Orbital_Charge& loc,
                                          const elecstate::ElecState* pelec,
                                          const psi::Psi<T>* psi,
                                          LCAO_Hamilt& uhm,
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
    //--------------------------------------------------------
    // implement four terms which needs integration
    //--------------------------------------------------------
    this->calForceStressIntegralPart(GlobalV::GAMMA_ONLY_LOCAL,
                                     isforce,
                                     isstress,
                                     loc,
                                     pelec,
                                     psi,
                                     foverlap,
                                     ftvnl_dphi,
                                     fvnl_dbeta,
                                     fvl_dphi,
                                     soverlap,
                                     stvnl_dphi,
                                     svnl_dbeta,
#ifdef __DEEPKS
                                     svl_dphi,
                                     svnl_dalpha,
#else
                                     svl_dphi,
#endif
                                     uhm,
                                     kv);
    // implement vdw force or stress here
    //  Peize Lin add 2014-04-04, update 2021-03-09
    //  jiyy add 2019-05-18, update 2021-05-02
    ModuleBase::matrix force_vdw;
    ModuleBase::matrix stress_vdw;
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
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

    // implement force from E-field
    ModuleBase::matrix fefield;
    if (GlobalV::EFIELD_FLAG && isforce)
    {
        fefield.create(nat, 3);
        elecstate::Efield::compute_force(GlobalC::ucell, fefield);
    }

    // implement force from E-field of tddft
    ModuleBase::matrix fefield_tddft;
    if (GlobalV::ESOLVER_TYPE == "TDDFT" && isforce)
    {
        fefield_tddft.create(nat, 3);
        elecstate::Efield::compute_force(GlobalC::ucell, fefield_tddft);
    }

    // implement force from gate field
    ModuleBase::matrix fgate;
    if (GlobalV::GATE_FLAG && isforce)
    {
        fgate.create(nat, 3);
        elecstate::Gatefield::compute_force(GlobalC::ucell, fgate);
    }
    // Force from implicit solvation model
    ModuleBase::matrix fsol;
    if (GlobalV::imp_sol && isforce)
    {
        fsol.create(nat, 3);
        GlobalC::solvent_model.cal_force_sol(GlobalC::ucell, rhopw, fsol);
    }
    // Force contribution from DFT+U
    ModuleBase::matrix force_dftu;
    ModuleBase::matrix stress_dftu;
    if (GlobalV::dft_plus_u) // Quxin add for DFT+U on 20201029
    {
        if (isforce)
        {
            force_dftu.create(nat, 3);
        }
        if (isstress)
        {
            stress_dftu.create(3, 3);
        }
        GlobalC::dftu.force_stress(pelec, *uhm.LM, force_dftu, stress_dftu, kv);
    }
    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->flk.finish_k();
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
                if (GlobalV::dft_plus_u)
                {
                    fcs(iat, i) += force_dftu(iat, i);
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
                if (GlobalV::EFIELD_FLAG)
                {
                    fcs(iat, i) += fefield(iat, i);
                }
                // E-field force of tddft
                if (GlobalV::ESOLVER_TYPE == "TDDFT")
                {
                    fcs(iat, i) += fefield_tddft(iat, i);
                }
                // Gate field force
                if (GlobalV::GATE_FLAG)
                {
                    fcs(iat, i) += fgate(iat, i);
                }
                // implicit solvation model
                if (GlobalV::imp_sol)
                {
                    fcs(iat, i) += fsol(iat, i);
                }
#ifdef __DEEPKS
                // mohan add 2021-08-04
                if (GlobalV::deepks_scf)
                {
                    fcs(iat, i) += GlobalC::ld.F_delta(iat, i);
                }
#endif
                // sum total force for correction
                sum += fcs(iat, i);
            }

            if (!(GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG))
            {
                for (int iat = 0; iat < nat; ++iat)
                {
                    fcs(iat, i) -= sum / nat;
                }
            }

            // xiaohui add "OUT_LEVEL", 2015-09-16
            if (GlobalV::OUT_LEVEL != "m")
            {
                GlobalV::ofs_running << " correction force for each atom along direction " << i + 1 << " is "
                                     << sum / nat << std::endl;
            }
        }

        if (GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG)
        {
            GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
        }

        // pengfei 2016-12-20
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            this->forceSymmetry(fcs, symm);
        }

#ifdef __DEEPKS
        // DeePKS force, caoyu add 2021-06-03
        if (GlobalV::deepks_out_labels) // not parallelized yet
        {
            GlobalC::ld.save_npy_f(fcs, "f_tot.npy", GlobalC::ucell.nat); // Ty/Bohr, F_tot
            if (GlobalV::deepks_scf)
            {
                GlobalC::ld.save_npy_f(fcs - GlobalC::ld.F_delta, "f_base.npy", GlobalC::ucell.nat); // Ry/Bohr, F_base

                if (GlobalV::GAMMA_ONLY_LOCAL)
                {
                    const std::vector<std::vector<double>>& dm_gamma
                        = dynamic_cast<const elecstate::ElecStateLCAO<double>*>(pelec)->get_DM()->get_DMK_vector();
                    GlobalC::ld.cal_gdmx(dm_gamma[0], GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, isstress);
                }
                else
                {
                    const std::vector<std::vector<std::complex<double>>>& dm_k
                        = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)
                              ->get_DM()
                              ->get_DMK_vector();
                    GlobalC::ld
                        .cal_gdmx_k(dm_k, GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, kv.nks, kv.kvec_d, isstress);
                }
                if (GlobalV::deepks_out_unittest)
                    GlobalC::ld.check_gdmx(GlobalC::ucell.nat);
                GlobalC::ld.cal_gvx(GlobalC::ucell.nat);

                if (GlobalV::deepks_out_unittest)
                    GlobalC::ld.check_gvx(GlobalC::ucell.nat);
                GlobalC::ld.save_npy_gvx(GlobalC::ucell.nat); //  /Bohr, grad_vx
            }
            else
            {
                GlobalC::ld.save_npy_f(fcs, "f_base.npy", GlobalC::ucell.nat); // no scf, F_base=F_tot
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
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "OVERLAP    FORCE", foverlap, 0);
            //  this->print_force("TVNL_DPHI  force",ftvnl_dphi,GlobalV::TEST_FORCE);
            //  this->print_force("VNL_DBETA  force",fvnl_dbeta,GlobalV::TEST_FORCE);
            // this->print_force("T_VNL      FORCE",ftvnl,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "T_VNL      FORCE", ftvnl, 0);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VL_dPHI    FORCE", fvl_dphi, 0);
            // this->print_force("VL_dPHI    FORCE",fvl_dphi,1,ry);
            // this->print_force("VL_dVL     FORCE",fvl_dvl,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VL_dVL     FORCE", fvl_dvl, 0);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EWALD      FORCE", fewalds, 0);
            // 	this->print_force("VLOCAL     FORCE",fvlocal,GlobalV::TEST_FORCE);
            // this->print_force("EWALD      FORCE",fewalds,1,ry);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "NLCC       FORCE", fcc, 0);
            ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "SCC        FORCE", fscc, 0);
            // this->print_force("NLCC       FORCE",fcc,1,ry);
            // this->print_force("SCC        FORCE",fscc,1,ry);
            //-------------------------------
            // put extra force here for test!
            //-------------------------------
            if (GlobalV::EFIELD_FLAG)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EFIELD     FORCE", fefield, 0);
                // this->print_force("EFIELD     FORCE",fefield,1,ry);
            }
            if (GlobalV::ESOLVER_TYPE == "TDDFT")
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "EFIELD_TDDFT     FORCE", fefield_tddft, 0);
                // this->print_force("EFIELD_TDDFT     FORCE",fefield_tddft,1,ry);
            }
            if (GlobalV::GATE_FLAG)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "GATEFIELD     FORCE", fgate, 0);
                // this->print_force("GATEFIELD     FORCE",fgate,1,ry);
            }
            if (GlobalV::imp_sol)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "IMP_SOL     FORCE", fsol, 0);
                // this->print_force("IMP_SOL     FORCE",fsol,1,ry);
            }
            if (vdw_solver != nullptr)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "VDW        FORCE", force_vdw, 0);
                // this->print_force("VDW        FORCE",force_vdw,1,ry);
            }
#ifdef __DEEPKS
            // caoyu add 2021-06-03
            if (GlobalV::deepks_scf)
            {
                ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "DeePKS 	FORCE", GlobalC::ld.F_delta, 1);
                // this->print_force("DeePKS 	FORCE", GlobalC::ld.F_delta, 1, ry);
            }
#endif
        }

        GlobalV::ofs_running << std::setiosflags(std::ios::left);

        // this->printforce_total(ry, istestf, fcs);
        ModuleIO::print_force(GlobalV::ofs_running, GlobalC::ucell, "TOTAL-FORCE (eV/Angstrom)", fcs, 0);
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
                if (GlobalV::dft_plus_u)
                {
                    scs(i, j) += stress_dftu(i, j);
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
            symm->symmetrize_mat3(scs, GlobalC::ucell);
        } // end symmetry

#ifdef __DEEPKS
        if (GlobalV::deepks_out_labels) // not parallelized yet
        {
            GlobalC::ld.save_npy_s(scs,
                                   "s_base.npy",
                                   GlobalC::ucell.omega); // change to energy unit Ry when printing, S_base;
        }
        if (GlobalV::deepks_scf)
        {
            if (ModuleSymmetry::Symmetry::symm_flag == 1)
            {
                symm->symmetrize_mat3(svnl_dalpha, GlobalC::ucell);
            } // end symmetry
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    scs(i, j) += svnl_dalpha(i, j);
                }
            }
        }
        if (GlobalV::deepks_out_labels) // not parallelized yet
        {
            // wenfei add 2021/11/2
            if (GlobalV::deepks_scf)
            {
                GlobalC::ld.save_npy_s(scs,
                                       "s_tot.npy",
                                       GlobalC::ucell.omega); // change to energy unit Ry when printing, S_tot, w/ model
                GlobalC::ld.cal_gvepsl(GlobalC::ucell.nat);
                GlobalC::ld.save_npy_gvepsl(GlobalC::ucell.nat); //  unitless, grad_vepsl
            }
            else
            {
                GlobalC::ld.save_npy_s(
                    scs,
                    "s_tot.npy",
                    GlobalC::ucell.omega); // change to energy unit Ry when printing, S_tot, w/o model;
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
            ModuleIO::print_stress("OVERLAP  STRESS", soverlap, GlobalV::TEST_STRESS, ry);
            // test
            ModuleIO::print_stress("T        STRESS", stvnl_dphi, GlobalV::TEST_STRESS, ry);
            ModuleIO::print_stress("VNL      STRESS", svnl_dbeta, GlobalV::TEST_STRESS, ry);

            ModuleIO::print_stress("T_VNL    STRESS", stvnl, GlobalV::TEST_STRESS, ry);

            ModuleIO::print_stress("VL_dPHI  STRESS", svl_dphi, GlobalV::TEST_STRESS, ry);
            ModuleIO::print_stress("VL_dVL   STRESS", sigmadvl, GlobalV::TEST_STRESS, ry);
            ModuleIO::print_stress("HAR      STRESS", sigmahar, GlobalV::TEST_STRESS, ry);

            ModuleIO::print_stress("EWALD    STRESS", sigmaewa, GlobalV::TEST_STRESS, ry);
            ModuleIO::print_stress("cc       STRESS", sigmacc, GlobalV::TEST_STRESS, ry);
            //		ModuleIO::print_stress("NLCC       STRESS",sigmacc,GlobalV::TEST_STRESS,ry);
            ModuleIO::print_stress("XC       STRESS", sigmaxc, GlobalV::TEST_STRESS, ry);
            if (vdw_solver != nullptr)
            {
                ModuleIO::print_stress("VDW      STRESS", sigmaxc, GlobalV::TEST_STRESS, ry);
            }
            if (GlobalV::dft_plus_u)
            {
                ModuleIO::print_stress("DFTU     STRESS", stress_dftu, GlobalV::TEST_STRESS, ry);
            }
            ModuleIO::print_stress("TOTAL    STRESS", scs, GlobalV::TEST_STRESS, ry);

        } // end of test
        GlobalV::ofs_running << std::setiosflags(std::ios::left);
        // print total stress
        ModuleIO::print_stress("TOTAL-STRESS", scs, true, ry);

        double unit_transform = 0.0;
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
        double external_stress[3] = {GlobalV::PRESS1, GlobalV::PRESS2, GlobalV::PRESS3};

        for (int i = 0; i < 3; i++)
        {
            scs(i, i) -= external_stress[i] / unit_transform;
        }
        GlobalV::PRESSURE = (scs(0, 0) + scs(1, 1) + scs(2, 2)) / 3;
    } // end of stress calculation

    ModuleBase::timer::tick("Force_Stress_LCAO", "getForceStress");
    return;
}

// local pseudopotential, ewald, core correction, scc terms in force
template<typename T>
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
    f_pw.cal_force_cc(fcc, rhopw, chr);
    //--------------------------------------------------------
    // force due to self-consistent charge.
    //--------------------------------------------------------
    f_pw.cal_force_scc(fscc, rhopw, vnew, vnew_exist);
    return;
}

// overlap, kinetic, nonlocal pseudopotential, Local potential terms in force and stress
template<>
void Force_Stress_LCAO<double>::calForceStressIntegralPart(const bool isGammaOnly,
    const bool isforce,
    const bool isstress,
    Local_Orbital_Charge& loc,
    const elecstate::ElecState* pelec,
    const psi::Psi<double>* psi,
    ModuleBase::matrix& foverlap,
    ModuleBase::matrix& ftvnl_dphi,
    ModuleBase::matrix& fvnl_dbeta,
    ModuleBase::matrix& fvl_dphi,
    ModuleBase::matrix& soverlap,
    ModuleBase::matrix& stvnl_dphi,
    ModuleBase::matrix& svnl_dbeta,
#if __DEEPKS
    ModuleBase::matrix& svl_dphi,
    ModuleBase::matrix& svnl_dalpha,
#else
    ModuleBase::matrix& svl_dphi,
#endif
    LCAO_Hamilt& uhm,
    const K_Vectors& kv)
{
    flk.ftable_gamma(isforce,
        isstress,
        psi,
        loc,
        pelec,
        foverlap,
        ftvnl_dphi,
        fvnl_dbeta,
        fvl_dphi,
        soverlap,
        stvnl_dphi,
        svnl_dbeta,
#if __DEEPKS
        svl_dphi,
        svnl_dalpha,
#else
        svl_dphi,
#endif
        uhm);
    return;
}
template<>
void Force_Stress_LCAO<std::complex<double>>::calForceStressIntegralPart(const bool isGammaOnly,
    const bool isforce,
    const bool isstress,
    Local_Orbital_Charge& loc,
    const elecstate::ElecState* pelec,
    const psi::Psi<std::complex<double>>* psi,
    ModuleBase::matrix& foverlap,
    ModuleBase::matrix& ftvnl_dphi,
    ModuleBase::matrix& fvnl_dbeta,
    ModuleBase::matrix& fvl_dphi,
    ModuleBase::matrix& soverlap,
    ModuleBase::matrix& stvnl_dphi,
    ModuleBase::matrix& svnl_dbeta,
#if __DEEPKS
    ModuleBase::matrix& svl_dphi,
    ModuleBase::matrix& svnl_dalpha,
#else
    ModuleBase::matrix& svl_dphi,
#endif
    LCAO_Hamilt& uhm,
    const K_Vectors& kv)
{
        flk.ftable_k(isforce,
                     isstress,
                     *this->RA,
                     psi,
                     loc,
                     pelec,
                     foverlap,
                     ftvnl_dphi,
                     fvnl_dbeta,
                     fvl_dphi,
                     soverlap,
                     stvnl_dphi,
                     svnl_dbeta,
#if __DEEPKS
                     svl_dphi,
                     svnl_dalpha,
#else
                     svl_dphi,
#endif
                     uhm,
            kv);
    return;
}

// vlocal, hartree, ewald, core correction, exchange-correlation terms in stress
template<typename T>
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
template<typename T>
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