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

void print_rhofft(ModulePW::PW_Basis* pw_rhod,
                  ModulePW::PW_Basis* pw_rho,
                  ModulePW::PW_Basis_Big* pw_big,
                  std::ofstream& ofs)
{
    std::cout << " UNIFORM GRID DIM        : " << pw_rho->nx << " * " << pw_rho->ny << " * " << pw_rho->nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG)   : " << pw_big->nbx << " * " << pw_big->nby << " * " << pw_big->nbz
              << std::endl;
    if (PARAM.globalv.double_grid)
    {
        std::cout << " UNIFORM GRID DIM(DENSE) : " << pw_rhod->nx << " * " << pw_rhod->ny << " * " << pw_rhod->nz
                  << std::endl;
    }

    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Setup plane waves of charge/potential:                          "
           "   |"
        << std::endl;
    ofs << " | Use the energy cutoff and the lattice vectors to generate the   "
           "   |"
        << std::endl;
    ofs << " | dimensions of FFT grid. The number of FFT grid on each "
           "processor   |"
        << std::endl;
    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space "
           "is   |"
        << std::endl;
    ofs << " | different for charege/potential and wave functions. We also set "
           "   |"
        << std::endl;
    ofs << " | the 'sticks' for the parallel of FFT. The number of plane waves "
           "   |"
        << std::endl;
    ofs << " | is 'npw' in each processor.                                     "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";
    ofs << "\n SETUP THE PLANE WAVE BASIS" << std::endl;

    double ecut = 4 * PARAM.inp.ecutwfc;
    if (PARAM.inp.nx * PARAM.inp.ny * PARAM.inp.nz > 0)
    {
        ecut = pw_rho->gridecut_lat * pw_rho->tpiba2;
        ofs << "use input fft dimensions for wave functions." << std::endl;
        ofs << "calculate energy cutoff from nx, ny, nz:" << std::endl;
    }

    ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for charge/potential (unit:Ry)", ecut);

    ModuleBase::GlobalFunc::OUT(ofs, "fft grid for charge/potential", pw_rho->nx, pw_rho->ny, pw_rho->nz);
    ModuleBase::GlobalFunc::OUT(ofs, "fft grid division", pw_big->bx, pw_big->by, pw_big->bz);
    ModuleBase::GlobalFunc::OUT(ofs, "big fft grid for charge/potential", pw_big->nbx, pw_big->nby, pw_big->nbz);
    ModuleBase::GlobalFunc::OUT(ofs, "nbxx", pw_big->nbxx);
    ModuleBase::GlobalFunc::OUT(ofs, "nrxx", pw_rho->nrxx);

    ofs << "\n SETUP PLANE WAVES FOR CHARGE/POTENTIAL" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", pw_rho->npwtot);
    ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", pw_rho->nstot);

    ofs << "\n PARALLEL PW FOR CHARGE/POTENTIAL" << std::endl;
    ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

    for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
    {
        ofs << " " << std::setw(8) << i + 1 << std::setw(15) << pw_rho->nst_per[i] << std::setw(15)
            << pw_rho->npw_per[i] << std::endl;
    }
    ofs << " --------------- sum -------------------" << std::endl;
    ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << pw_rho->nstot << std::setw(15)
        << pw_rho->npwtot << std::endl;

    ModuleBase::GlobalFunc::OUT(ofs, "number of |g|", pw_rho->ngg);
    ModuleBase::GlobalFunc::OUT(ofs, "max |g|", pw_rho->gg_uniq[pw_rho->ngg - 1]);
    ModuleBase::GlobalFunc::OUT(ofs, "min |g|", pw_rho->gg_uniq[0]);

    if (PARAM.globalv.double_grid)
    {
        ofs << std::endl;
        ofs << std::endl;
        ofs << std::endl;
        double ecut = PARAM.inp.ecutrho;
        if (PARAM.inp.ndx * PARAM.inp.ndy * PARAM.inp.ndz > 0)
        {
            ecut = pw_rhod->gridecut_lat * pw_rhod->tpiba2;
            ofs << "use input fft dimensions for the dense part of charge "
                   "density."
                << std::endl;
            ofs << "calculate energy cutoff from ndx, ndy, ndz:" << std::endl;
        }
        ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for dense charge/potential (unit:Ry)", ecut);

        ModuleBase::GlobalFunc::OUT(ofs, "fft grid for dense charge/potential", pw_rhod->nx, pw_rhod->ny, pw_rhod->nz);

        ModuleBase::GlobalFunc::OUT(ofs, "nrxx", pw_rhod->nrxx);

        ofs << "\n SETUP PLANE WAVES FOR dense CHARGE/POTENTIAL" << std::endl;
        ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", pw_rhod->npwtot);
        ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", pw_rhod->nstot);

        ofs << "\n PARALLEL PW FOR dense CHARGE/POTENTIAL" << std::endl;
        ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

        for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
        {
            ofs << " " << std::setw(8) << i + 1 << std::setw(15) << pw_rhod->nst_per[i] << std::setw(15)
                << pw_rhod->npw_per[i] << std::endl;
        }
        ofs << " --------------- sum -------------------" << std::endl;
        ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << pw_rhod->nstot << std::setw(15)
            << pw_rhod->npwtot << std::endl;

        ModuleBase::GlobalFunc::OUT(ofs, "number of |g|", pw_rhod->ngg);
        ModuleBase::GlobalFunc::OUT(ofs, "max |g|", pw_rhod->gg_uniq[pw_rhod->ngg - 1]);
        ModuleBase::GlobalFunc::OUT(ofs, "min |g|", pw_rhod->gg_uniq[0]);
    }
}

void print_wfcfft(const Input_para& inp, ModulePW::PW_Basis_K& pw_wfc, std::ofstream& ofs)
{
    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Setup plane waves of wave functions:                            "
           "   |"
        << std::endl;
    ofs << " | Use the energy cutoff and the lattice vectors to generate the   "
           "   |"
        << std::endl;
    ofs << " | dimensions of FFT grid. The number of FFT grid on each "
           "processor   |"
        << std::endl;
    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space "
           "is   |"
        << std::endl;
    ofs << " | different for charege/potential and wave functions. We also set "
           "   |"
        << std::endl;
    ofs << " | the 'sticks' for the parallel of FFT. The number of plane wave "
           "of  |"
        << std::endl;
    ofs << " | each k-point is 'npwk[ik]' in each processor                    "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";
    ofs << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;

    double ecut = inp.ecutwfc;
    if (std::abs(ecut - pw_wfc.gk_ecut * pw_wfc.tpiba2) > 1e-6)
    {
        ecut = pw_wfc.gk_ecut * pw_wfc.tpiba2;
        ofs << "Energy cutoff for wavefunc is incompatible with nx, ny, nz and "
               "it will be reduced!"
            << std::endl;
    }
    ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for wavefunc (unit:Ry)", ecut);
    ModuleBase::GlobalFunc::OUT(ofs, "fft grid for wave functions", pw_wfc.nx, pw_wfc.ny, pw_wfc.nz);
    ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", pw_wfc.npwtot);
    ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", pw_wfc.nstot);

    ofs << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
    ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

    for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
    {
        ofs << " " << std::setw(8) << i + 1 << std::setw(15) << pw_wfc.nst_per[i] << std::setw(15) << pw_wfc.npw_per[i]
            << std::endl;
    }

    ofs << " --------------- sum -------------------" << std::endl;
    ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << pw_wfc.nstot << std::setw(15)
        << pw_wfc.npwtot << std::endl;
    ModuleBase::GlobalFunc::DONE(ofs, "INIT PLANEWAVE");
}

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