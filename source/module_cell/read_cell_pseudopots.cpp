#include <cstring> // Peize Lin fix bug about strcmp 2016-08-02

#include "module_base/parallel_common.h"
#include "module_io/input.h"
#include "unitcell.h"

//==========================================================
// Read pseudopotential according to the dir
//==========================================================
void UnitCell::read_cell_pseudopots(const std::string &pp_dir, std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell","read_cell_pseudopots");
	// setup reading log for pseudopot_upf
	std::stringstream ss;
	ss << GlobalV::global_out_dir << "atom_pseudo.log";
	
	// Read in the atomic pseudo potentials
	std::string pp_address;
	for (int i = 0;i < ntype;i++)
	{
		Pseudopot_upf upf;
	
		// mohan update 2010-09-12	
		int error = 0;
		int error_ap = 0;
		
		if(GlobalV::MY_RANK==0)
		{
			pp_address = pp_dir + this->pseudo_fn[i];
			error = upf.init_pseudo_reader( pp_address, this->pseudo_type[i] ); //xiaohui add 2013-06-23

			if(error==0) // mohan add 2021-04-16
			{
				if(this->atoms[i].flag_empty_element)	// Peize Lin add for bsse 2021.04.07
				{
					upf.set_empty_element();			
				}
                upf.set_upf_q(); // liuyu add 2023-09-21
                upf.check_atwfc_norm(); // liuyu add 2023-09-26
                // average pseudopotential if needed
                error_ap = upf.average_p(GlobalV::soc_lambda); // added by zhengdy 2020-10-20
            }
        }

#ifdef __MPI
		Parallel_Common::bcast_int(error);
		Parallel_Common::bcast_int(error_ap);
#endif

		if(error_ap) 
		{
            ModuleBase::WARNING_QUIT("UnitCell::read_cell_pseudopots", "error when average the pseudopotential.");
        }

        if (error == 1)
        {
            std::cout << " Pseudopotential directory now is : " << pp_address << std::endl;
            GlobalV::ofs_warning << " Pseudopotential directory now is : " << pp_address << std::endl;
            ModuleBase::WARNING_QUIT("UnitCell::read_cell_pseudopots", "Couldn't find pseudopotential file.");
        }
        else if (error == 2)
        {
            ModuleBase::WARNING_QUIT("UnitCell::read_cell_pseudopots", "Pseudopotential data do not match.");
        }
        else if (error == 3)
        {
            ModuleBase::WARNING_QUIT(
                "UnitCell::read_cell_pseudopots",
                "Check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave "
                "functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n "
                "unreasonable large (should be near 1.0), ABACUS would quit. \n The solution is to turn off the wave "
                "functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
        }
        else if (error == 4)
        {
            ModuleBase::WARNING_QUIT("UnitCell::read_cell_pseudopots", "Unknown pseudopotential type.");
        }

        if (GlobalV::MY_RANK == 0)
        {
            atoms[i].ncpp.set_pseudo(upf);

            log << "\n Read in pseudopotential file is " << pseudo_fn[i] << std::endl;
            ModuleBase::GlobalFunc::OUT(log, "pseudopotential type", atoms[i].ncpp.pp_type);
            ModuleBase::GlobalFunc::OUT(log, "exchange-correlation functional", atoms[i].ncpp.xc_func);
            ModuleBase::GlobalFunc::OUT(log, "nonlocal core correction", atoms[i].ncpp.nlcc);
            // ModuleBase::GlobalFunc::OUT(log, "spin orbital", atoms[i].has_so);
            ModuleBase::GlobalFunc::OUT(log, "valence electrons", atoms[i].ncpp.zv);
            ModuleBase::GlobalFunc::OUT(log, "lmax", atoms[i].ncpp.lmax);
            ModuleBase::GlobalFunc::OUT(log, "number of zeta", atoms[i].ncpp.nchi);
            ModuleBase::GlobalFunc::OUT(log, "number of projectors", atoms[i].ncpp.nbeta);
            for (int ib = 0; ib < atoms[i].ncpp.nbeta; ib++)
            {
                ModuleBase::GlobalFunc::OUT(log, "L of projector", atoms[i].ncpp.lll[ib]);
            }
            //			ModuleBase::GlobalFunc::OUT(log,"Grid Mesh Number", atoms[i].mesh);
            if (GlobalV::DFT_FUNCTIONAL != "default")
            {
                std::string xc_func1 = GlobalV::DFT_FUNCTIONAL;
                transform(xc_func1.begin(), xc_func1.end(), xc_func1.begin(), (::toupper));
                if (xc_func1 != atoms[i].ncpp.xc_func)
                {
                    std::cout << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << std::endl;
                    std::cout << " dft_functional in pseudopot file is: " << atoms[i].ncpp.xc_func << std::endl;
                    std::cout << " Please make sure this is what you need" << std::endl;
                    GlobalV::ofs_warning << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << std::endl;
                    GlobalV::ofs_warning << " dft_functional in pseudopot file is: " << atoms[i].ncpp.xc_func
                                         << std::endl;
                    GlobalV::ofs_warning << " Please make sure this is what you need" << std::endl;

                    atoms[i].ncpp.xc_func = xc_func1;
                    log << " XC functional updated to : " << GlobalV::DFT_FUNCTIONAL << std::endl;
                    ModuleBase::GlobalFunc::OUT(log, "exchange-correlation functional", atoms[i].ncpp.xc_func);
                }
            }
        }
    }
    return;
}


void UnitCell::print_unitcell_pseudo(const std::string &fn)
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("UnitCell","print_unitcell_pseudo");
	std::ofstream ofs( fn.c_str() );

	this->print_cell(ofs);
	for (int i = 0;i < ntype;i++)
	{
		atoms[i].print_Atom(ofs);
	}

	ofs.close();
	return;
}
