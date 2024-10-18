#include "print_info.h"

#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"


namespace ModuleIO
{

void setup_parameters(UnitCell& ucell, K_Vectors& kv)
{
    ModuleBase::TITLE("ModuleIO", "setup_parameters");

    if(PARAM.inp.calculation=="scf" || PARAM.inp.calculation=="relax" || PARAM.inp.calculation=="cell-relax" || PARAM.inp.calculation=="nscf"
	        || PARAM.inp.calculation=="get_pchg" || PARAM.inp.calculation=="get_wf" || PARAM.inp.calculation=="md")
	{
		std::cout << " ---------------------------------------------------------" << std::endl;
		if(PARAM.inp.calculation=="scf")
		{
			std::cout << " Self-consistent calculations for electrons" << std::endl;
		}
		else if(PARAM.inp.calculation=="test")
		{
			std::cout << " Test run" << std::endl;
		}
		if(PARAM.inp.calculation=="relax")
		{
            std::cout << " Ion relaxation calculations" << std::endl;
		}
        if(PARAM.inp.calculation=="cell-relax")
        {
            std::cout << " Cell relaxation calculations" << std::endl;
        }
		if(PARAM.inp.calculation=="md")
		{
			std::cout << " Molecular Dynamics simulations" << std::endl;

			std::cout << " ---------------------------------------------------------" << std::endl;

            if (PARAM.mdp.md_type == "fire")
            {
                std::cout << " ENSEMBLE                 : " << "FIRE" << std::endl;
            }
            else if (PARAM.mdp.md_type == "nve")
            {
                std::cout << " ENSEMBLE                 : " << "NVE" << std::endl;
            }
            else if (PARAM.mdp.md_type == "nvt")
            {
                std::cout << " ENSEMBLE                 : "
                          << "NVT    mode: " << PARAM.mdp.md_thermostat << std::endl;
            }
            else if (PARAM.mdp.md_type == "npt")
            {
                std::cout << " ENSEMBLE                 : "
                          << "NPT    mode: " << PARAM.mdp.md_pmode << std::endl;
            }
            else if (PARAM.mdp.md_type == "langevin")
            {
                std::cout << " ENSEMBLE                 : " << "Langevin" << std::endl;
            }
            else if (PARAM.mdp.md_type == "msst")
            {
                std::cout << " ENSEMBLE                 : " << "MSST" << std::endl;
            }

            std::cout << " Time interval(fs)        : " << PARAM.mdp.md_dt << std::endl;
        }
        std::cout << " ---------------------------------------------------------" << std::endl;


		std::cout << " " << std::setw(8) << "SPIN"
		     << std::setw(16) << "KPOINTS"
		     << std::setw(12) << "PROCESSORS"
             << std::setw(12) << "THREADS";

		const bool orbinfo = (PARAM.inp.basis_type=="lcao" || PARAM.inp.basis_type=="lcao_in_pw" 
						  || (PARAM.inp.basis_type=="pw" && PARAM.inp.init_wfc.substr(0, 3) == "nao"));
		if (orbinfo) { std::cout << std::setw(12) << "NBASE"; }

		std::cout << std::endl;
		std::cout << " " << std::setw(8) << PARAM.inp.nspin;

		if(PARAM.globalv.gamma_only_local)
		{
			std::cout << std::setw(16) << "Gamma";
		}
		else
		{
			std::cout << std::setw(16) << kv.get_nkstot();
		}

		std::cout << std::setw(12) << GlobalV::NPROC
		     << std::setw(12) << PARAM.globalv.nthread_per_proc * GlobalV::NPROC;
		if (orbinfo) { std::cout << std::setw(12) << PARAM.globalv.nlocal; }

		std::cout << std::endl;




		std::cout << " ---------------------------------------------------------" << std::endl;
		if(PARAM.inp.basis_type == "lcao")
		{
			std::cout << " Use Systematically Improvable Atomic bases" << std::endl;
		}
		else if(PARAM.inp.basis_type == "lcao_in_pw")
		{
			std::cout << " Expand Atomic bases into plane waves" << std::endl;
		}
		else if(PARAM.inp.basis_type == "pw")
		{
			std::cout << " Use plane wave basis" << std::endl;
		}
		std::cout << " ---------------------------------------------------------" << std::endl;



		//----------------------------------
		// second part
		//----------------------------------

		std::cout << " " << std::setw(8) << "ELEMENT";

		if (orbinfo)
		{
			std::cout << std::setw(16) << "ORBITALS";
			std::cout << std::setw(12) << "NBASE";
		}
		std::cout << std::setw(12) << "NATOM";

		std::cout << std::setw(12) << "XC";
		std::cout << std::endl;


		const std::string spectrum = "spdfghi";
		for(int it=0; it<ucell.ntype; ++it)
		{
			std::cout << " " << std::setw(8) << ucell.atoms[it].label;

			if (orbinfo)
			{
				std::stringstream orb;
				int norb = 0;

				for(int L=0; L<=ucell.atoms[it].nwl; ++L)        // pengfei Li 16-2-29
				{
					norb += (2*L+1)* ucell.atoms[it].l_nchi[L];
					orb << ucell.atoms[it].l_nchi[L];
					orb << spectrum[L];
				}
				orb << "-" << ucell.atoms[it].Rcut << "au";
				
				std::cout << std::setw(16) << orb.str();
				std::cout << std::setw(12) << norb;
			}


			std::cout << std::setw(12) << ucell.atoms[it].na;
			std::cout << std::endl;
		}

		std::cout << " ---------------------------------------------------------" << std::endl;
		std::cout << " Initial plane wave basis and FFT box" << std::endl;
		std::cout << " ---------------------------------------------------------" << std::endl;

	}

	return;
}

void print_time(time_t& time_start, time_t& time_finish)
{
    // print out information before ABACUS ends
	std::cout << "\n START  Time  : " << ctime(&time_start);
	std::cout << " FINISH Time  : " << ctime(&time_finish);
	std::cout << " TOTAL  Time  : " << int(difftime(time_finish, time_start)) << std::endl;
	std::cout << " SEE INFORMATION IN : " << PARAM.globalv.global_out_dir << std::endl;

	GlobalV::ofs_running << "\n Start  Time  : " << ctime(&time_start);
	GlobalV::ofs_running << " Finish Time  : " << ctime(&time_finish);

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	GlobalV::ofs_running << " Total  Time  : " << unsigned(hour) << " h "
	    << unsigned(mins) << " mins "
	    << unsigned(secs) << " secs "<< std::endl;
}

/*
void ModuleIO::print_scf(const int &istep, const int &iter)
{
    if(PARAM.inp.basis_type=="pw")
    {
        GlobalV::ofs_running << "\n PW ALGORITHM ------------- ";
    }
    else
    {
        GlobalV::ofs_running << "\n LCAO ALGORITHM ------------- ";
    }

    if(PARAM.inp.calculation=="scf")
    {
        GlobalV::ofs_running << "ELEC = " << std::setw(4) << unsigned(iter);
    }
    else if(PARAM.inp.calculation=="relax" || PARAM.inp.calculation=="cell-relax")
    {
        GlobalV::ofs_running << "ION = " << std::setw(4) << unsigned(istep+1)
                             << "  ELEC = " << std::setw(4) << unsigned(iter);
    }
    else if(PARAM.inp.calculation=="md")
    {
        GlobalV::ofs_running << "MD = " << std::setw(4) << unsigned(istep+1)
                             << "  ELEC = " << std::setw(4) << unsigned(iter);
    }

    GlobalV::ofs_running << " --------------------------------\n";
}
*/

void print_screen(const int& stress_step, const int& force_step, const int& istep)
{
    std::cout << " -------------------------------------------" << std::endl;
	GlobalV::ofs_running << "\n -------------------------------------------" << std::endl;

	if(PARAM.inp.calculation=="scf") //add 4 lines 2015-09-06, xiaohui
	{
        std::cout << " SELF-CONSISTENT : " << std::endl;
		GlobalV::ofs_running << " SELF-CONSISTENT" << std::endl;
	}
	else if(PARAM.inp.calculation=="nscf") //add 4 lines 2015-09-06, xiaohui
	{
        std::cout << " NONSELF-CONSISTENT : " << std::endl;
		GlobalV::ofs_running << " NONSELF-CONSISTENT" << std::endl;
	}
	else if(PARAM.inp.calculation=="md")
	{
        std::cout << " STEP OF MOLECULAR DYNAMICS : " << unsigned(istep) << std::endl;
		GlobalV::ofs_running << " STEP OF MOLECULAR DYNAMICS : " << unsigned(istep) << std::endl;
	}
	else
	{
		if(PARAM.inp.relax_new)
		{
			std::cout << " STEP OF RELAXATION : " << unsigned(istep) << std::endl;
			GlobalV::ofs_running << " STEP OF RELAXATION : " << unsigned(istep) << std::endl;
		}
		else if(PARAM.inp.calculation=="relax") //pengfei 2014-10-13
		{
        	std::cout << " STEP OF ION RELAXATION : " << unsigned(istep) << std::endl;
			GlobalV::ofs_running << " STEP OF ION RELAXATION : " << unsigned(istep) << std::endl;
		}
    	else if(PARAM.inp.calculation=="cell-relax")
    	{
        	std::cout << " RELAX CELL : " << unsigned(stress_step) << std::endl;
        	std::cout << " RELAX IONS : " << unsigned(force_step) << " (in total: " << unsigned(istep) << ")" << std::endl;
			GlobalV::ofs_running << " RELAX CELL : " << unsigned(stress_step) << std::endl;
        	GlobalV::ofs_running << " RELAX IONS : " << unsigned(force_step) << " (in total: " << unsigned(istep) << ")" << std::endl;
    	}
	}

    std::cout << " -------------------------------------------" << std::endl;
    GlobalV::ofs_running << " -------------------------------------------" << std::endl;
}

} // namespace ModuleIO