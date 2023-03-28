#include "relax_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h" // use chr.
#include "module_io/print_info.h"
#include "module_io/write_wfc_r.h"

void Relax_Driver::relax_driver(ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("Ions","opt_ions");
	ModuleBase::timer::tick("Ions","opt_ions");

	if(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
	{
		if(!GlobalV::relax_new)
		{
			rl_old.init_relax();
		}
		else
		{
			rl.init_relax(GlobalC::ucell.nat, Lattice_Change_Basic::out_stru);
		}
	}

    this->istep = 1;
	int force_step = 1;           // pengfei Li 2018-05-14
	int stress_step = 1;
	bool stop= false;
	
    while(istep <= GlobalV::RELAX_NMAX && !stop)
    {
		time_t estart = time(NULL);

		if(GlobalV::OUT_LEVEL=="ie" &&
			(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax"
			|| GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="nscf"))
		{
			Print_Info::print_screen(stress_step, force_step, istep);
		}

		// mohan added eiter to count for the electron iteration number, 2021-01-28	
		p_esolver->Run(istep-1,GlobalC::ucell);

		time_t eend = time(NULL);
		time_t fstart = time(NULL);

        if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
			//I'm considering putting force and stress
			//as part of ucell and use ucell to pass information
			//back and forth between esolver and relaxation
			//but I'll use force and stress explicitly here for now

            // calculate the total energy
            p_esolver->cal_Energy(GlobalC::en.etot);

			//calculate and gather all parts of total ionic forces
			ModuleBase::matrix force;
			if(GlobalV::CAL_FORCE)
			{
				p_esolver->cal_Force(force);
			}
			//calculate and gather all parts of stress
			ModuleBase::matrix stress;
			if(GlobalV::CAL_STRESS)
			{
				p_esolver->cal_Stress(stress);
			}

			if(GlobalV::relax_new)
			{
				stop = rl.relax_step(force, stress, GlobalC::en.etot);
			}
			else
			{
				stop = rl_old.relax_step(force, stress, istep, force_step, stress_step);    // pengfei Li 2018-05-14
			}
			
		}
		time_t fend = time(NULL);

		++istep;

    }

	if(GlobalV::OUT_LEVEL=="i")
	{
		std::cout << " ION DYNAMICS FINISHED :)" << std::endl;
	}

	ModuleBase::timer::tick("Ions","opt_ions_pw");
    return;
}