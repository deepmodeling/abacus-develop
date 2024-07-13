#include "relax_driver.h"

#include "module_base/global_file.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h" // use chr.
#include "module_io/json_output/output_info.h"
#include "module_io/print_info.h"
#include "module_io/write_wfc_r.h"

void Relax_Driver::relax_driver(ModuleESolver::ESolver* p_esolver)
{
    ModuleBase::TITLE("Ions", "opt_ions");
    ModuleBase::timer::tick("Ions", "opt_ions");

    if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
    {
        if (!GlobalV::relax_new)
        {
            rl_old.init_relax(GlobalC::ucell.nat);
        }
        else
        {
            rl.init_relax(GlobalC::ucell.nat);
        }
    }

    this->istep = 1;
    int force_step = 1; // pengfei Li 2018-05-14
    int stress_step = 1;
    bool stop = false;

    while (istep <= GlobalV::RELAX_NMAX && !stop)
    {
        time_t estart = time(nullptr);

        if (GlobalV::OUT_LEVEL == "ie"
            && (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax" || GlobalV::CALCULATION == "scf"
                || GlobalV::CALCULATION == "nscf")
            && (GlobalV::ESOLVER_TYPE != "lr"))
        {
            Print_Info::print_screen(stress_step, force_step, istep);
        }

#ifdef __RAPIDJSON
        Json::init_output_array_obj();
#endif //__RAPIDJSON

        // mohan added eiter to count for the electron iteration number, 2021-01-28
        p_esolver->runner(istep - 1, GlobalC::ucell);

        time_t eend = time(nullptr);
        time_t fstart = time(nullptr);
        ModuleBase::matrix force;
        ModuleBase::matrix stress;
        if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            // I'm considering putting force and stress
            // as part of ucell and use ucell to pass information
            // back and forth between esolver and relaxation
            // but I'll use force and stress explicitly here for now

            // calculate the total energy
            this->etot = p_esolver->cal_energy();

            // calculate and gather all parts of total ionic forces
            if (GlobalV::CAL_FORCE)
            {
                p_esolver->cal_force(force);
            }
            // calculate and gather all parts of stress
            if (GlobalV::CAL_STRESS)
            {
                p_esolver->cal_stress(stress);
            }

            if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
            {
                if (GlobalV::relax_new)
                {
                    stop = rl.relax_step(force, stress, this->etot);
                }
                else
                {
                    stop = rl_old.relax_step(istep,
                                             this->etot,
                                             GlobalC::ucell,
                                             force,
                                             stress,
                                             force_step,
                                             stress_step); // pengfei Li 2018-05-14
                }
                // print structure
                // changelog 20240509
                // because I move out the dependence on GlobalV from UnitCell::print_stru_file
                // so its parameter is calculated here
                bool need_orb = GlobalV::BASIS_TYPE == "pw";
                need_orb = need_orb && GlobalV::psi_initializer;
                need_orb = need_orb && GlobalV::init_wfc.substr(0, 3) == "nao";
                need_orb = need_orb || GlobalV::BASIS_TYPE == "lcao";
                need_orb = need_orb || GlobalV::BASIS_TYPE == "lcao_in_pw";
                std::stringstream ss, ss1;
                ss << GlobalV::global_out_dir << "STRU_ION_D";
                GlobalC::ucell.print_stru_file(ss.str(),
                                               GlobalV::NSPIN,
                                               true,
                                               GlobalV::CALCULATION == "md",
                                               GlobalV::out_mul,
                                               need_orb,
                                               GlobalV::deepks_setorb,
                                               GlobalV::MY_RANK);

                if (Ions_Move_Basic::out_stru)
                {
                    ss1 << GlobalV::global_out_dir << "STRU_ION";
                    ss1 << istep << "_D";
                    GlobalC::ucell.print_stru_file(ss1.str(),
                                                   GlobalV::NSPIN,
                                                   true,
                                                   GlobalV::CALCULATION == "md",
                                                   GlobalV::out_mul,
                                                   need_orb,
                                                   GlobalV::deepks_setorb,
                                                   GlobalV::MY_RANK);

                    GlobalC::ucell.print_cell_cif("STRU_NOW.cif");
                }

                if (p_esolver && stop && p_esolver->get_maxniter() == p_esolver->get_niter()
                    && !(p_esolver->get_conv_elec()))
                {
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                              << std::endl;
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                              << std::endl;
                    std::cout << " Relaxation is converged, but the SCF is unconverged! The results are unreliable. "
                              << std::endl;
                    std::cout
                        << " It is suggested to increase the maximum SCF step and/or perform the relaxation again."
                        << std::endl;
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                              << std::endl;
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                              << std::endl;
                    GlobalV::ofs_running
                        << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running
                        << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running
                        << "\n Relaxation is converged, but the SCF is unconverged! The results are unreliable.. "
                        << std::endl;
                    GlobalV::ofs_running
                        << "\n It is suggested to increase the maximum SCF step and/or perform the relaxation again. "
                        << std::endl;
                    GlobalV::ofs_running
                        << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running
                        << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                }
            }

#ifdef __RAPIDJSON
            // add the energy to outout
            Json::add_output_energy(p_esolver->cal_energy() * ModuleBase::Ry_to_eV);
#endif
        }
#ifdef __RAPIDJSON
        // add Json of cell coo stress force
        double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
        double fac = ModuleBase::Ry_to_eV / 0.529177;
        Json::add_output_cell_coo_stress_force(&GlobalC::ucell, force, fac, stress, unit_transform);
#endif //__RAPIDJSON

        time_t fend = time(nullptr);

        ++istep;
    }

    if (GlobalV::OUT_LEVEL == "i")
    {
        std::cout << " ION DYNAMICS FINISHED :)" << std::endl;
    }

    if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
    {
        ModuleBase::Global_File::delete_tmp_files();
    }

    ModuleBase::timer::tick("Ions", "opt_ions");
    return;
}
