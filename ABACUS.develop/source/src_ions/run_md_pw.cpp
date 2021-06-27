#include "src_ions/run_md_pw.h"
#include "../src_pw/forces.h"
#include "../src_pw/stress_pw.h"
#include "../src_pw/global.h" // use chr.
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../src_pw/pw_complement.h"
#include "../src_pw/pw_basis.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "src_ions/MD_basic.h"

void Run_MD_PW::md_ions_pw(void)
{
    TITLE("Run_MD_PW", "md_ions_pw");
    timer::tick("Run_MD_PW", "md_ions_pw", 'C');

    if (OUT_LEVEL == "i")
    {
        cout << setprecision(12);
        cout << " " << setw(7) << "ISTEP"
             << setw(5) << "NE"
             << setw(15) << "ETOT(eV)"
             << setw(15) << "EDIFF(eV)"
             << setw(15) << "MAX_F(eV/A)"
             << setw(15) << "TRADIUS(Bohr)"
             << setw(8) << "UPDATE"
             << setw(11) << "ETIME(MIN)"
             << setw(11) << "FTIME(MIN)"
             << endl;
    }

    // allocation for ion movement.
    if (FORCE)
    {
        IMM.allocate();
        CE.allocate_ions();
    }

    if (STRESS) // pengfei Li 2018-05-14
    {
        LCM.allocate();
    }

    MD_basic mdb(INPUT.mdp, ucell);
    int mdtype = INPUT.mdp.mdtype;

    this->istep = 1;
    int force_step = 1; // pengfei Li 2018-05-14
    int stress_step = 1;
    bool stop = false;

    while (istep <= NSTEP && !stop)
    {
        time_t estart = time(NULL);

        if (OUT_LEVEL == "ie")
        {
            cout << " -------------------------------------------" << endl;    
            cout << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
            cout << " -------------------------------------------" << endl;
            ofs_running << " -------------------------------------------" << endl;
            ofs_running << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
            ofs_running << " -------------------------------------------" << endl;
        }

        if (vdwd2_para.flag_vdwd2) //Peize Lin add 2014-04-03, update 2021-03-09
        {
            Vdwd2 vdwd2(ucell, vdwd2_para);
            vdwd2.cal_energy();
            en.evdw = vdwd2.get_energy();
        }
        if (vdwd3_para.flag_vdwd3) //jiyy add 2019-05-18, update 2021-05-02
        {
            Vdwd3 vdwd3(ucell, vdwd3_para);
            vdwd3.cal_energy();
            en.evdw = vdwd3.get_energy();
        }

        // mohan added eiter to count for the electron iteration number, 2021-01-28
        int eiter = 0;
        if (CALCULATION == "md")
        {
            if (Exx_Global::Hybrid_Type::No == exx_global.info.hybrid_type)
            {
                elec.self_consistent(istep - 1);
                eiter = elec.iter;
            }
            else if (Exx_Global::Hybrid_Type::Generate_Matrix == exx_global.info.hybrid_type)
            {
                throw invalid_argument(TO_STRING(__FILE__) + TO_STRING(__LINE__));
            }
            else // Peize Lin add 2019-03-09
            {
                if (exx_global.info.separate_loop)
                {
                    for (size_t hybrid_step = 0; hybrid_step != exx_global.info.hybrid_step; ++hybrid_step)
                    {
                        elec.self_consistent(istep - 1);
                        eiter += elec.iter;
                        if (elec.iter == 1 || hybrid_step == exx_global.info.hybrid_step - 1) // exx converge
                            break;
                        exx_global.info.set_xcfunc(xcf);
                        exx_lip.cal_exx();
                    }
                }
                else
                {
                    elec.self_consistent(istep - 1);
                    eiter += elec.iter;
                    exx_global.info.set_xcfunc(xcf);
                    elec.self_consistent(istep - 1);
                    eiter += elec.iter;
                }
            }
        }
        // mohan added 2021-01-28, perform stochastic calculations
        else if (CALCULATION == "md-sto")
        {
            elec_sto.scf_stochastic(istep - 1);
            eiter = elec_sto.iter;
        }

        CE.update_all_pos(ucell);

        if (pot.out_potential == 2)
        {
            stringstream ssp;
            stringstream ssp_ave;
            ssp << global_out_dir << "ElecStaticPot";
            ssp_ave << global_out_dir << "ElecStaticPot_AVE";
            pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        time_t eend = time(NULL);
        time_t fstart = time(NULL);

        if (mdtype == 1 || mdtype == 2)
        {
            mdb.runNVT(istep);
        }
        else if (mdtype == 0)
        {
            mdb.runNVE(istep);
        }
        else if (mdtype == -1)
        {
            stop = mdb.runFIRE(istep);
        }
        else
        {
            WARNING_QUIT("opt_ions", "mdtype should be -1~2!");
        }

        time_t fend = time(NULL);

        CE.save_pos_next(ucell);

        //xiaohui add CE.istep = istep 2014-07-07
        CE.update_istep(istep);

        // charge extrapolation if istep>0.
        CE.extrapolate_charge();

        //reset local potential and initial wave function
        pot.init_pot(istep, pw.strucFac);
        ofs_running << " Setup the new wave functions?" << endl;
        wf.wfcinit();

        if (OUT_LEVEL == "i")
        {
            double etime_min = difftime(eend, estart) / 60.0;
            double ftime_min = difftime(fend, fstart) / 60.0;
            stringstream ss;
            ss << MOVE_IONS << istep;

            cout << " " << setw(7) << ss.str()
                 << setw(5) << eiter
                 << setw(15) << setprecision(6) << en.etot * Ry_to_eV
                 << setw(15) << IMM.get_ediff() * Ry_to_eV
                 << setprecision(3)
                 << setw(15) << IMM.get_largest_grad() * Ry_to_eV / 0.529177
                 << setw(15) << IMM.get_trust_radius()
                 << setw(8) << IMM.get_update_iter()
                 << setprecision(2) << setw(11) << etime_min
                 << setw(11) << ftime_min << endl;
        }

        ++istep;
    }

    if (CALCULATION == "scf" || CALCULATION == "relax" || CALCULATION == "cell-relax")
    {
        ofs_running << "\n\n --------------------------------------------" << endl;
        ofs_running << setprecision(16);
        ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl;
        ofs_running << " --------------------------------------------\n\n"
                    << endl;
    }

    if (OUT_LEVEL == "i")
    {
        cout << " ION DYNAMICS FINISHED :)" << endl;
    }

    timer::tick("Run_MD_PW", "md_ions_pw", 'C');
    return;
}

void Run_MD_PW::md_cells_pw()
{
    TITLE("Run_MD_PW", "md_cells_pw");
    timer::tick("Run_MD_PW", "md_cells_pw", 'C');

    wf.allocate(kv.nks);

    UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    ppcell.init(ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    hm.hpw.allocate(wf.npwx, NPOL, ppcell.nkb, pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc(pw.nggm, ppcell.vloc);
    DONE(ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    ppcell.init_vnl(ucell);
    DONE(ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    pot.init_pot(0, pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    pot.newd();

    DONE(ofs_running, "INIT POTENTIAL");

    //==================================================
    // create ppcell.tab_at , for trial wave functions.
    //==================================================
    wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if (NBANDS != 0 || (CALCULATION != "scf-sto" && CALCULATION != "relax-sto" && CALCULATION != "md-sto")) //qianrui add
    {
        wf.wfcinit();
    }

    switch (exx_global.info.hybrid_type) // Peize Lin add 2019-03-09
    {
    case Exx_Global::Hybrid_Type::HF:
    case Exx_Global::Hybrid_Type::PBE0:
    case Exx_Global::Hybrid_Type::HSE:
        exx_lip.init(&kv, &wf, &pw, &UFFT, &ucell);
        break;
    case Exx_Global::Hybrid_Type::No:
        break;
    case Exx_Global::Hybrid_Type::Generate_Matrix:
    default:
        throw invalid_argument(TO_STRING(__FILE__) + TO_STRING(__LINE__));
    }

    DONE(ofs_running, "INIT BASIS");

    // ion optimization begins
    // electron density optimization is included in ion optimization

    this->md_ions_pw();

    ofs_running << "\n\n --------------------------------------------" << endl;
    ofs_running << setprecision(16);
    ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl;
    ofs_running << " --------------------------------------------\n\n" << endl;

    timer::tick("Run_MD_PW", "md_cells_pw", 'C');
}