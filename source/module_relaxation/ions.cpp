#include "ions.h"
#include "../src_pw/global.h" // use chr.
#include "../src_io/print_info.h"
#include "variable_cell.h" // mohan add 2021-02-01
#include "src_io/write_wfc_realspace.h"

void Ions::opt_ions(ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("Ions","opt_ions");
	ModuleBase::timer::tick("Ions","opt_ions");
	
	if(GlobalV::OUT_LEVEL=="i")
	{
		std::cout << std::setprecision(12);
    	std::cout<< " " << std::setw(7)<< "ISTEP" 
		<<std::setw(5)<< "NE"
		<<std::setw(15)<< "ETOT(eV)"
		<<std::setw(15)<< "EDIFF(eV)"
        <<std::setw(15)<< "MAX_F(eV/A)"
        <<std::setw(15)<< "TRADIUS(Bohr)"
		<<std::setw(8)<< "UPDATE"
		<<std::setw(11)<< "ETIME(MIN)"
		<<std::setw(11)<< "FTIME(MIN)"
        <<std::endl;
	}

    // Geometry optimization algorithm setup.
    if (GlobalV::CALCULATION=="relax")
    {
        //Ions_Move_Methods
        IMM.allocate();
    }
    if (GlobalV::CALCULATION=="cell-relax")
    {
        //Ions_Move_Methods
        IMM.allocate();
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    } 

    this->istep = 1;
	int force_step = 1;           // pengfei Li 2018-05-14
	int stress_step = 1;
	bool stop= false;
	
    while(istep <= GlobalV::RELAX_NMAX && !stop)
    {
		time_t estart = time(NULL);

		if(GlobalV::OUT_LEVEL=="ie")
		{
			Print_Info::print_screen(stress_step, force_step, istep);
		}

		// mohan added eiter to count for the electron iteration number, 2021-01-28	
		p_esolver->Run(istep-1,GlobalC::ucell);

		time_t eend = time(NULL);
		time_t fstart = time(NULL);

        if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax" || GlobalV::CALCULATION.substr(0,3)=="sto")
        {
			//I'm considering putting force and stress
			//as part of ucell and use ucell to pass information
			//back and forth between esolver and relaxation
			//but I'll use force and stress explicitly here for now
			
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
			stop = this->relaxation(force, stress, istep, force_step, stress_step);    // pengfei Li 2018-05-14
		}
		time_t fend = time(NULL);

		if(GlobalV::OUT_LEVEL=="i")
		{
			double etime_min = difftime(eend, estart)/60.0; 
			double ftime_min = difftime(fend, fstart)/60.0; 
			std::stringstream ss;
			ss << GlobalV::RELAX_METHOD << istep;
			
			std::cout << " " << std::setw(7) << ss.str() 
			<< std::setw(5) << p_esolver->getniter()
			<< std::setw(15) << std::setprecision(6) << GlobalC::en.etot * ModuleBase::Ry_to_eV 
			<< std::setw(15) << IMM.get_ediff() * ModuleBase::Ry_to_eV
			<< std::setprecision(3)
			<< std::setw(15) << IMM.get_largest_grad() * ModuleBase::Ry_to_eV / 0.529177
			<< std::setw(15) << IMM.get_trust_radius()
			<< std::setw(8) << IMM.get_update_iter()
			<< std::setprecision(2) << std::setw(11) << etime_min
			<< std::setw(11) << ftime_min << std::endl;
		}

		++istep;

    }

	if(GlobalV::OUT_LEVEL=="i")
	{
		std::cout << " ION DYNAMICS FINISHED :)" << std::endl;
	}

	ModuleBase::timer::tick("Ions","opt_ions_pw");
    return;
}