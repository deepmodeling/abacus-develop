#include "FORCE_STRESS_tddft.h"

#include "../src_pw/global.h"
#include "./dftu.h" //Quxin add for DFT+U on 20201029
// new
#include "../module_base/timer.h"
#include "../module_surchem/efield.h" // liuyu add 2022-05-18
#include "../module_surchem/gatefield.h" // liuyu add 2022-09-13
#include "../module_surchem/surchem.h" //sunml add 2022-08-10
#include "../src_pw/forces.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h" //caoyu add for deepks 2021-06-03
#endif

double Force_Stress_LCAO_TDDFT::force_invalid_threshold_ev = 0.00;
double Force_Stress_LCAO_TDDFT::output_acc = 1.0e-8;

Force_Stress_LCAO_TDDFT::Force_Stress_LCAO_TDDFT(Record_adj& ra) : RA(&ra)
{
}
Force_Stress_LCAO_TDDFT::~Force_Stress_LCAO_TDDFT()
{
}

void Force_Stress_LCAO_TDDFT::getForceStress(const bool isforce,
                                             const bool isstress,
                                             const bool istestf,
                                             const bool istests,
                                             Local_Orbital_Charge& loc,
                                             const psi::Psi<double>* psid,
                                             const psi::Psi<std::complex<double>>* psi,
                                             LCAO_Hamilt& uhm,
                                             ModuleBase::matrix& fcs,
                                             ModuleBase::matrix& scs,
                                             ModuleBase::Vector3<double>* vel)
{
    ModuleBase::TITLE("Force_Stress_LCAO_TDDFT", "getForceStress");
    ModuleBase::timer::tick("Force_Stress_LCAO_TDDFT", "getForceStress");

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
    ModuleBase::matrix ftddft;

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
        ftddft.create(nat, 3);
        // calculate basic terms in Force, same method with PW base
        this->calForcePwPart(fvl_dvl, fewalds, fcc, fscc);
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
        this->calStressPwPart(sigmadvl, sigmahar, sigmaewa, sigmacc, sigmaxc);
    }
    //--------------------------------------------------------
    // implement four terms which needs integration
    //--------------------------------------------------------
    this->calForceStressIntegralPart(GlobalV::GAMMA_ONLY_LOCAL,
                                     isforce,
                                     isstress,
                                     loc,
                                     psid,
                                     psi,
                                     foverlap,
                                     ftddft,
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
                                     vel);
    // implement vdw force or stress here
    //  Peize Lin add 2014-04-04, update 2021-03-09
    ModuleBase::matrix force_vdw;
    ModuleBase::matrix stress_vdw;
    if (GlobalC::vdwd2_para.flag_vdwd2)
    {
        if (isforce)
        {
            force_vdw.create(nat, 3);
            Vdwd2 vdwd2(GlobalC::ucell, GlobalC::vdwd2_para);
            vdwd2.cal_force();
            for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
            {
                force_vdw(iat, 0) = vdwd2.get_force()[iat].x;
                force_vdw(iat, 1) = vdwd2.get_force()[iat].y;
                force_vdw(iat, 2) = vdwd2.get_force()[iat].z;
            }
        }
        if (isstress)
        {
            Vdwd2 vdwd2(GlobalC::ucell, GlobalC::vdwd2_para);
            vdwd2.cal_stress();
            stress_vdw = vdwd2.get_stress().to_matrix();
        }
    }
    // jiyy add 2019-05-18, update 2021-05-02
    else if (GlobalC::vdwd3_para.flag_vdwd3)
    {
        if (isforce)
        {
            force_vdw.create(nat, 3);
            Vdwd3 vdwd3(GlobalC::ucell, GlobalC::vdwd3_para);
            vdwd3.cal_force();
            for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
            {
                force_vdw(iat, 0) = vdwd3.get_force()[iat].x;
                force_vdw(iat, 1) = vdwd3.get_force()[iat].y;
                force_vdw(iat, 2) = vdwd3.get_force()[iat].z;
            }
        }
        if (isstress)
        {
            Vdwd3 vdwd3(GlobalC::ucell, GlobalC::vdwd3_para);
            vdwd3.cal_stress();
            stress_vdw = vdwd3.get_stress().to_matrix();
        }
    }
    // implement force from E-field
    ModuleBase::matrix fefield;
    if (GlobalV::EFIELD_FLAG && isforce)
    {
        fefield.create(nat, 3);
        Efield::compute_force(GlobalC::ucell, fefield);
    }
    // implement force from gate field
    ModuleBase::matrix fgate;
    if (GlobalV::GATE_FLAG && isforce)
    {
        fgate.create(nat, 3);
        Gatefield::compute_force(GlobalC::ucell, fgate);
    }
    // Force from implicit solvation model
    ModuleBase::matrix fsol;
    if (GlobalV::imp_sol && isforce)
    {
        fsol.create(nat, 3);
        GlobalC::solvent_model.cal_force_sol(GlobalC::ucell, GlobalC::rhopw, fsol);
    }
    // Force contribution from DFT+U
    ModuleBase::matrix force_dftu;
    ModuleBase::matrix stress_dftu;
    if (INPUT.dft_plus_u)
    {
        // Quxin add for DFT+U on 20201029
        GlobalC::dftu.force_stress(loc.dm_gamma, loc.dm_k, *uhm.LM);

        if (isforce)
        {
            force_dftu.create(nat, 3);
        }
        if (isstress)
        {
            stress_dftu.create(3, 3);
        }
        for (int i = 0; i < 3; i++)
        {
            if (isstress)
            {
                for (int j = 0; j < 3; j++)
                {
                    stress_dftu(j, i) = GlobalC::dftu.stress_dftu.at(j).at(i);
                }
            }
            if (isforce)
            {
                for (int iat = 0; iat < nat; iat++)
                {
                    force_dftu(iat, i) = GlobalC::dftu.force_dftu.at(iat).at(i);
                }
            }
        }
    }
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
                //fcs(iat, i)+= foverlap(iat, i) + ftvnl_dphi(iat, i)+ fvnl_dbeta(iat, i)
                 fcs(iat, i) += foverlap(iat, i) + ftddft(iat, i) + ftvnl_dphi(iat, i) + fvnl_dbeta(iat, i)
                       + fvl_dphi(iat, i) + fvl_dvl(iat, i) // derivative of local potential force (pw)
                       + fewalds(iat, i) // ewald force (pw)
                       + fcc(iat, i) // nonlinear core correction force (pw)
                       + fscc(iat, i); // self consistent corretion force (pw)

                // Force contribution from DFT+U, Quxin add on 20201029
                if (INPUT.dft_plus_u)
                {
                    fcs(iat, i) += force_dftu(iat, i);
                }
                // VDW force of vdwd2 or vdwd3
                if (GlobalC::vdwd2_para.flag_vdwd2 || GlobalC::vdwd3_para.flag_vdwd3)
                {
                    fcs(iat, i) += force_vdw(iat, i);
                }
                // E-field force
                if (GlobalV::EFIELD_FLAG)
                {
                    fcs(iat, i) += fefield(iat, i);
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
        if (ModuleSymmetry::Symmetry::symm_flag)
        {
            this->forceSymmetry(fcs);
        }

#ifdef __DEEPKS
        // DeePKS force, caoyu add 2021-06-03
        if (GlobalV::deepks_out_labels) // not parallelized yet
        {
            const Parallel_Orbitals* pv = loc.ParaV;
            GlobalC::ld.save_npy_f(fcs, "f_tot.npy", GlobalC::ucell.nat); // Ty/Bohr, F_tot
            if (GlobalV::deepks_scf)
            {
                GlobalC::ld.save_npy_f(fcs - GlobalC::ld.F_delta, "f_base.npy", GlobalC::ucell.nat); // Ry/Bohr, F_base

                if (GlobalV::GAMMA_ONLY_LOCAL)
                {
                    GlobalC::ld.cal_gdmx(loc.dm_gamma[0],
                                         GlobalC::ucell,
                                         GlobalC::ORB,
                                         GlobalC::GridD,
                                         pv->trace_loc_row,
                                         pv->trace_loc_col,
                                         isstress);
                }
                else
                {
                    GlobalC::ld.cal_gdmx_k(loc.dm_k,
                                           GlobalC::ucell,
                                           GlobalC::ORB,
                                           GlobalC::GridD,
                                           pv->trace_loc_row,
                                           pv->trace_loc_col,
                                           GlobalC::kv.nks,
                                           GlobalC::kv.kvec_d,
                                           isstress);
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
            GlobalV::ofs_running << std::setiosflags(ios::showpos);
            GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(8) << std::endl;
            //-----------------------------
            // regular force terms test.
            //-----------------------------
            // this->print_force("OVERLAP    FORCE",foverlap,1,ry);
            f_pw.print("OVERLAP    FORCE", foverlap, 0);
            f_pw.print("TDDFT    FORCE", ftddft, 0);

            //  this->print_force("TVNL_DPHI  force",ftvnl_dphi,GlobalV::TEST_FORCE);
            //  this->print_force("VNL_DBETA  force",fvnl_dbeta,GlobalV::TEST_FORCE);
            // this->print_force("T_VNL      FORCE",ftvnl,1,ry);
            f_pw.print("T_VNL      FORCE", ftvnl, 0);
            f_pw.print("VL_dPHI    FORCE", fvl_dphi, 0);
            // this->print_force("VL_dPHI    FORCE",fvl_dphi,1,ry);
            // this->print_force("VL_dVL     FORCE",fvl_dvl,1,ry);
            f_pw.print("VL_dVL     FORCE", fvl_dvl, 0);
            f_pw.print("EWALD      FORCE", fewalds, 0);
            // 	this->print_force("VLOCAL     FORCE",fvlocal,GlobalV::TEST_FORCE);
            // this->print_force("EWALD      FORCE",fewalds,1,ry);
            f_pw.print("NLCC       FORCE", fcc, 0);
            f_pw.print("SCC        FORCE", fscc, 0);
            // this->print_force("NLCC       FORCE",fcc,1,ry);
            // this->print_force("SCC        FORCE",fscc,1,ry);
            //-------------------------------
            // put extra force here for test!
            //-------------------------------
            if (GlobalV::EFIELD_FLAG)
            {
                f_pw.print("EFIELD     FORCE", fefield, 0);
                // this->print_force("EFIELD     FORCE",fefield,1,ry);
            }
            if (GlobalV::GATE_FLAG)
            {
                f_pw.print("GATEFIELD     FORCE", fgate, 0);
                // this->print_force("GATEFIELD     FORCE",fgate,1,ry);
            }
            if (GlobalV::imp_sol)
            {
                f_pw.print("IMP_SOL     FORCE", fsol, 0);
                // this->print_force("IMP_SOL     FORCE",fsol,1,ry);
            }
            if (GlobalC::vdwd2_para.flag_vdwd2 || GlobalC::vdwd3_para.flag_vdwd3)
            {
                f_pw.print("VDW        FORCE", force_vdw, 0);
                // this->print_force("VDW        FORCE",force_vdw,1,ry);
            }
#ifdef __DEEPKS
            // caoyu add 2021-06-03
            if (GlobalV::deepks_scf)
            {
                this->print_force("DeePKS 	FORCE", GlobalC::ld.F_delta, 1, ry);
            }
#endif
        }

        GlobalV::ofs_running << std::setiosflags(ios::left);

        // this->printforce_total(ry, istestf, fcs);
        f_pw.print("   TOTAL-FORCE (eV/Angstrom)", fcs, 0);
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
                    if (abs(fcs(iat, i) * ModuleBase::Ry_to_eV / 0.529177)
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
                             + sigmadvl(i, j) // derivative of local potential stress (pw)
                             + sigmaewa(i, j) // ewald stress (pw)
                             + sigmacc(i, j) // nonlinear core correction stress (pw)
                             + sigmaxc(i, j) // exchange corretion stress
                             + sigmahar(i, j); // hartree stress

                // VDW stress from linpz and jiyy
                if (GlobalC::vdwd2_para.flag_vdwd2 || GlobalC::vdwd3_para.flag_vdwd3)
                {
                    scs(i, j) += stress_vdw(i, j);
                }
                // DFT plus U stress from qux
                if (INPUT.dft_plus_u)
                {
                    scs(i, j) += stress_dftu(i, j);
                }
            }
        }

#ifdef __DEEPKS
        if (GlobalV::deepks_out_labels) // not parallelized yet
        {
            GlobalC::ld.save_npy_s(scs,
                                   "s_base.npy",
                                   GlobalC::ucell.omega); // change to energy unit Ry when printing, S_base;
            // wenfei add 2021/11/2
            if (GlobalV::deepks_scf)
            {
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        scs(i, j) += svnl_dalpha(i, j);
                    }
                }
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

        if (ModuleSymmetry::Symmetry::symm_flag)
        {
            GlobalC::symm.stress_symmetry(scs, GlobalC::ucell);
        } // end symmetry

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
            GlobalV::ofs_running << std::setiosflags(ios::showpos);
            GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(10) << std::endl;
            sc_pw.print_stress("OVERLAP  STRESS", soverlap, GlobalV::TEST_STRESS, ry);
            // test
            sc_pw.print_stress("T        STRESS", stvnl_dphi, GlobalV::TEST_STRESS, ry);
            sc_pw.print_stress("VNL      STRESS", svnl_dbeta, GlobalV::TEST_STRESS, ry);

            sc_pw.print_stress("T_VNL    STRESS", stvnl, GlobalV::TEST_STRESS, ry);

            sc_pw.print_stress("VL_dPHI  STRESS", svl_dphi, GlobalV::TEST_STRESS, ry);
            sc_pw.print_stress("VL_dVL   STRESS", sigmadvl, GlobalV::TEST_STRESS, ry);
            sc_pw.print_stress("HAR      STRESS", sigmahar, GlobalV::TEST_STRESS, ry);

            sc_pw.print_stress("EWALD    STRESS", sigmaewa, GlobalV::TEST_STRESS, ry);
            sc_pw.print_stress("cc       STRESS", sigmacc, GlobalV::TEST_STRESS, ry);
            //		sc_pw.print_stress("NLCC       STRESS",sigmacc,GlobalV::TEST_STRESS,ry);
            sc_pw.print_stress("XC       STRESS", sigmaxc, GlobalV::TEST_STRESS, ry);
            if (GlobalC::vdwd2_para.flag_vdwd2 || GlobalC::vdwd3_para.flag_vdwd3)
            {
                sc_pw.print_stress("VDW      STRESS", sigmaxc, GlobalV::TEST_STRESS, ry);
            }
            if (INPUT.dft_plus_u)
            {
                sc_pw.print_stress("DFTU     STRESS", sigmaxc, GlobalV::TEST_STRESS, ry);
            }
            sc_pw.print_stress("TOTAL    STRESS", scs, GlobalV::TEST_STRESS, ry);

        } // end of test
        GlobalV::ofs_running << std::setiosflags(ios::left);
        // print total stress
        sc_pw.printstress_total(scs, ry);

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

// overlap, kinetic, nonlocal pseudopotential, Local potential terms in force and stress
void Force_Stress_LCAO_TDDFT::calForceStressIntegralPart(const bool isGammaOnly,
                                                         const bool isforce,
                                                         const bool isstress,
                                                         Local_Orbital_Charge& loc,
                                                         const psi::Psi<double>* psid,
                                                         const psi::Psi<std::complex<double>>* psi,
                                                         ModuleBase::matrix& foverlap,
                                                         ModuleBase::matrix& ftddft,
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
                                                         ModuleBase::Vector3<double>* vel)
{
    if (isGammaOnly)
    {
        flk.ftable_gamma(isforce,
                         isstress,
                         psid,
                         loc,
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
    }
    else
    {
        flk_tddft.ftable_k(isforce,
                           isstress,
                           *this->RA,
                           psi,
                           loc,
                           foverlap,
                           ftddft,
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
                           vel);
    }
    return;
}
