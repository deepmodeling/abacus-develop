#include "relax_driver.h"

#include "module_base/global_file.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h" // use chr.
#include "module_io/cif_io.h"
#include "module_io/json_output/output_info.h"
#include "module_io/output_log.h"
#include "module_io/print_info.h"
#include "module_io/read_exit_file.h"
#include "module_io/write_wfc_r.h"
#include "module_parameter/parameter.h"
void Relax_Driver::relax_driver(ModuleESolver::ESolver* p_esolver, UnitCell& ucell)
{
    ModuleBase::TITLE("Ions", "opt_ions");
    ModuleBase::timer::tick("Ions", "opt_ions");

    if (PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax" )
    {
        if (!PARAM.inp.relax_new)
        {
            rl_old.init_relax(ucell.nat);
        }
        else
        {
            rl.init_relax(ucell.nat);
        }
    }
    
    this->istep = 1;
    int force_step = 1; // pengfei Li 2018-05-14
    int stress_step = 1;
    bool stop = false;

    while (istep <= PARAM.inp.relax_nmax && !stop)
    {
        time_t estart = time(nullptr);

        if (PARAM.inp.out_level == "ie"
            && (PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax" || PARAM.inp.calculation == "scf"
                || PARAM.inp.calculation == "nscf")
            && (PARAM.inp.esolver_type != "lr"))
        {
            ModuleIO::print_screen(stress_step, force_step, istep);
        }

#ifdef __RAPIDJSON
        Json::init_output_array_obj();
#endif //__RAPIDJSON

        // mohan added eiter to count for the electron iteration number, 2021-01-28
        p_esolver->runner(ucell, istep - 1);

        time_t eend = time(nullptr);
        time_t fstart = time(nullptr);
        ModuleBase::matrix force;
        ModuleBase::matrix stress;
        if (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax")
        {
            // I'm considering putting force and stress
            // as part of ucell and use ucell to pass information
            // back and forth between esolver and relaxation
            // but I'll use force and stress explicitly here for now

            // calculate the total energy
            this->etot = p_esolver->cal_energy();

            // calculate and gather all parts of total ionic forces
            if (PARAM.inp.cal_force)
            {
                p_esolver->cal_force(ucell, force);
            }
            // calculate and gather all parts of stress
            if (PARAM.inp.cal_stress)
            {
                p_esolver->cal_stress(ucell, stress);
            }

            if (PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax")
            {
                if (PARAM.inp.relax_new)
                {
                    stop = rl.relax_step(ucell,force, stress, this->etot);
                }
                else
                {
                    stop = rl_old.relax_step(istep,
                                             this->etot,
                                             ucell,
                                             force,
                                             stress,
                                             force_step,
                                             stress_step); // pengfei Li 2018-05-14
                }         
                // print structure
                // changelog 20240509
                // because I move out the dependence on GlobalV from UnitCell::print_stru_file
                // so its parameter is calculated here
                bool need_orb = PARAM.inp.basis_type == "pw";
                need_orb = need_orb && PARAM.inp.psi_initializer;
                need_orb = need_orb && PARAM.inp.init_wfc.substr(0, 3) == "nao";
                need_orb = need_orb || PARAM.inp.basis_type == "lcao";
                need_orb = need_orb || PARAM.inp.basis_type == "lcao_in_pw";
                std::stringstream ss, ss1;
                ss << PARAM.globalv.global_out_dir << "STRU_ION_D";
                ucell.print_stru_file(ss.str(),
                                      PARAM.inp.nspin,
                                      true,
                                      PARAM.inp.calculation == "md",
                                      PARAM.inp.out_mul,
                                      need_orb,
                                      PARAM.globalv.deepks_setorb,
                                      GlobalV::MY_RANK);

                if (Ions_Move_Basic::out_stru)
                {
                    ss1 << PARAM.globalv.global_out_dir << "STRU_ION";
                    ss1 << istep << "_D";
                    ucell.print_stru_file(ss1.str(),
                                          PARAM.inp.nspin,
                                          true,
                                          PARAM.inp.calculation == "md",
                                          PARAM.inp.out_mul,
                                          need_orb,
                                          PARAM.globalv.deepks_setorb,
                                          GlobalV::MY_RANK);
                    ModuleIO::CifParser::write(PARAM.globalv.global_out_dir + "STRU_NOW.cif",
                                               ucell,
                                               "# Generated by ABACUS ModuleIO::CifParser",
                                               "data_?");
                }

                ModuleIO::output_after_relax(stop, p_esolver->conv_esolver, GlobalV::ofs_running);
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
        Json::add_output_cell_coo_stress_force(&ucell, force, fac, stress, unit_transform);
#endif //__RAPIDJSON

        if (stop == false)
        {
            stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, "EXIT", GlobalV::ofs_running);
        }

        time_t fend = time(nullptr);

        ++istep;
    }

    if (PARAM.inp.out_level == "i")
    {
        std::cout << " ION DYNAMICS FINISHED :)" << std::endl;
    }

    ModuleBase::timer::tick("Ions", "opt_ions");
    return;
}
