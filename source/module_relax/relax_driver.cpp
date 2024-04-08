#include "relax_driver.h"

#include "module_base/global_file.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h" // use chr.
#include "module_io/print_info.h"
#include "module_io/write_wfc_r.h"

#include "module_io/json_output/output_info.h"



template<typename FPTYPE, typename Device>
void Relax_Driver<FPTYPE, Device>::relax_driver(ModuleESolver::ESolver *p_esolver)
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
        time_t estart = time(NULL);

        if (GlobalV::OUT_LEVEL == "ie"
            && (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax" || GlobalV::CALCULATION == "scf"
                || GlobalV::CALCULATION == "nscf"))
        {
            Print_Info::print_screen(stress_step, force_step, istep);
        }

#ifdef __RAPIDJSON
        Json::init_output_array_obj();
#endif //__RAPIDJSON 

        // mohan added eiter to count for the electron iteration number, 2021-01-28
        p_esolver->run(istep - 1, GlobalC::ucell);

        time_t eend = time(NULL);
        time_t fstart = time(NULL);
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
                std::stringstream ss, ss1;
                ss << GlobalV::global_out_dir << "STRU_ION_D";
                GlobalC::ucell.print_stru_file(ss.str(), 2, 0);

                if (Ions_Move_Basic::out_stru)
                {
                    ss1 << GlobalV::global_out_dir << "STRU_ION";
                    ss1 << istep << "_D";
                    GlobalC::ucell.print_stru_file(ss1.str(), 2, 0);

                    GlobalC::ucell.print_cell_cif("STRU_NOW.cif");
                }

                ModuleESolver::ESolver_KS<FPTYPE, Device>* p_esolver_ks 
                = dynamic_cast<ModuleESolver::ESolver_KS<FPTYPE, Device>*>(p_esolver);
                if (p_esolver_ks 
                    && stop 
                    && p_esolver_ks->maxniter == p_esolver_ks->niter 
                    && !(p_esolver_ks->conv_elec))
                {
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    std::cout << " Relaxation is converged, but the SCF is unconverged! The results are unreliable. " << std::endl;
                    std::cout << " It is suggested to increase the maximum SCF step and/or perform the relaxation again." << std::endl;
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running << "\n Relaxation is converged, but the SCF is unconverged! The results are unreliable.. " << std::endl;
                    GlobalV::ofs_running << "\n It is suggested to increase the maximum SCF step and/or perform the relaxation again. " << std::endl;
                    GlobalV::ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                    GlobalV::ofs_running << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                }
            }
        }
#ifdef __RAPIDJSON
        //add Json of cell coo stress force
        double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
        double fac = ModuleBase::Ry_to_eV / 0.529177;
        Json::add_output_cell_coo_stress_force(
            &GlobalC::ucell,
            force,fac,
            stress,unit_transform);
#endif //__RAPIDJSON 
    
        time_t fend = time(NULL);

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

template class Relax_Driver<float, psi::DEVICE_CPU>;
template class Relax_Driver<double, psi::DEVICE_CPU>;
