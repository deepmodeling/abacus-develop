#include "relax_driver.h"

#include "source_base/global_file.h"
#include "source_pw/module_pwdft/global.h" // use chr.
#include "source_io/cif_io.h"
#include "source_io/json_output/output_info.h"
#include "source_io/output_log.h"
#include "source_io/print_info.h"
#include "source_io/read_exit_file.h"
#include "source_io/module_parameter/parameter.h"
#include "source_cell/print_cell.h"

void Relax_Driver::relax_driver(
		ModuleESolver::ESolver* p_esolver, 
		UnitCell& ucell,
		const Input_para& inp)
{ 
    ModuleBase::TITLE("Relax_Driver", "relax_driver");
    ModuleBase::timer::tick("Relax_Driver", "relax_driver");

    if (inp.calculation == "relax" || inp.calculation == "cell-relax" )
    {
        if (!inp.relax_new) // traditional relax
        {
            rl_old.init_relax(ucell.nat);
        }
        else // relax new
        {
            rl.init_relax(ucell.nat);
        }
    }

    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;

    while (istep <= inp.relax_nmax && !stop)
    {
        time_t estart = time(nullptr);

		if (inp.out_level == "ie"
				&& (inp.calculation == "relax" 
					|| inp.calculation == "cell-relax" 
					|| inp.calculation == "scf"
					|| inp.calculation == "nscf")
				&& (inp.esolver_type != "lr"))
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

        // I'm considering putting force and stress
        // as part of ucell and use ucell to pass information
        // back and forth between esolver and relaxation
        // but I'll use force and stress explicitly here for now

        // calculate the total energy
        this->etot = p_esolver->cal_energy();

        // calculate and gather all parts of total ionic forces
        if (inp.cal_force)
        {
            p_esolver->cal_force(ucell, force);
        }
		else
		{
            // do nothing
		}


        // calculate and gather all parts of stress
        if (inp.cal_stress)
        {
            p_esolver->cal_stress(ucell, stress);
        }
		else
		{
            // do nothing
		}

        if (inp.calculation == "relax" || inp.calculation == "cell-relax")
        {
            if (inp.relax_new)
            {
                stop = rl.relax_step(ucell, force, stress, this->etot);
                // mohan added 2025-07-14
                stress_step = istep+1;
                force_step = 1;
            }
            else
            {
                stop = rl_old.relax_step(istep,
                                         this->etot,
                                         ucell,
                                         force,
                                         stress,
                                         force_step,
                                         stress_step);
            }

            bool need_orb = inp.basis_type == "pw";
            need_orb = need_orb && inp.init_wfc.substr(0, 3) == "nao";
            need_orb = need_orb || inp.basis_type == "lcao";
            need_orb = need_orb || inp.basis_type == "lcao_in_pw";

            std::stringstream ss, ss1;
            ss << PARAM.globalv.global_out_dir << "STRU_ION_D";

            unitcell::print_stru_file(ucell,
                                  ucell.atoms,
                                  ucell.latvec,
                                  ss.str(),
                                  inp.nspin,
                                  true,
                                  inp.calculation == "md",
                                  inp.out_mul,
                                  need_orb,
                                  PARAM.globalv.deepks_setorb,
                                  GlobalV::MY_RANK);

            if (Ions_Move_Basic::out_stru)
            {
                ss1 << PARAM.globalv.global_out_dir << "STRU_ION";
                ss1 << istep << "_D";
                unitcell::print_stru_file(ucell,
                                      ucell.atoms,
                                      ucell.latvec,
                                      ss1.str(),
                                      inp.nspin,
                                      true,
                                      inp.calculation == "md",
                                      inp.out_mul,
                                      need_orb,
                                      PARAM.globalv.deepks_setorb,
                                      GlobalV::MY_RANK);
                ModuleIO::CifParser::write(PARAM.globalv.global_out_dir + "STRU_NOW.cif",
                                           ucell,
                                           "# Generated by ABACUS ModuleIO::CifParser",
                                           "data_?");
            }

            ModuleIO::output_after_relax(stop, p_esolver->conv_esolver, GlobalV::ofs_running);
        }// end relax or cell_relax

#ifdef __RAPIDJSON
        // add the energy to outout
        Json::add_output_energy(p_esolver->cal_energy() * ModuleBase::Ry_to_eV);
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
    } // end while (istep <= inp.relax_nmax && !stop)


	if (inp.calculation == "relax" || inp.calculation == "cell-relax")
	{
		if (istep-1 == inp.relax_nmax)
		{
			std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
			std::cout << " Geometry relaxation stops here due to reaching the maximum      " << std::endl;
			std::cout << " relaxation steps. More steps are needed to converge the results " << std::endl;
			std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
		}
		else
		{
			std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
			std::cout << " Geometry relaxation thresholds are reached within " << istep-1 << " steps." << std::endl; 
			std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
		}
	}
	else
	{
		// do nothing
	}

	if (inp.relax_nmax == 0)
    {
        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << " relax_nmax = 0, DRY RUN TEST SUCCEEDS :)" << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;
    }
	else
	{
		// do nothing 
	}

    ModuleBase::timer::tick("Relax_Driver", "relax_driver");
    return;
}
