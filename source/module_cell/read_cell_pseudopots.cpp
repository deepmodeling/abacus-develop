#include "unitcell_pseudo.h"
#include "../src_parallel/parallel_common.h"
#include "../input.h"
#ifdef __LCAO
//#include "../module_orbital/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif


#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

//==========================================================
// Read pseudopotential according to the dir
//==========================================================
void UnitCell_pseudo::read_cell_pseudopots(const std::string &pp_dir, std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell_pseudo","read_cell_pseudopots");
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
			error = upf.init_pseudo_reader( pp_address ); //xiaohui add 2013-06-23

			if(error==0) // mohan add 2021-04-16
			{
				if(this->atoms[i].flag_empty_element)	// Peize Lin add for bsse 2021.04.07
				{
					upf.set_empty_element();			
				}
				//average pseudopotential if needed
				error_ap = upf.average_p(GlobalV::soc_lambda); //added by zhengdy 2020-10-20
			}
		}

#ifdef __MPI
		Parallel_Common::bcast_int(error);
		Parallel_Common::bcast_int(error_ap);
#endif

		if(error_ap) 
		{
			ModuleBase::WARNING_QUIT("UnitCell_pseudo::read_pseudopot","error when average the pseudopotential.");
		}

		if(error==1)
		{
			std::cout << " Pseudopotential directory now is : " << pp_address << std::endl;
			GlobalV::ofs_warning << " Pseudopotential directory now is : " << pp_address << std::endl;
			ModuleBase::WARNING_QUIT("read_pseudopot","Couldn't find pseudopotential file.");
		}
		else if(error==2)
		{
			ModuleBase::WARNING_QUIT("read_pseudopot","Pseudopotential data do not match.");
		}
		else if(error==3)
		{
			ModuleBase::WARNING_QUIT("read_pseudopot","Check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n unreasonable large (should be near 1.0), ABACUS would quit. \n The solution is to turn off the wave functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
		}

//xiaohui add 2015-03-24
#ifdef __MPI
		Parallel_Common::bcast_bool(upf.functional_error);
#endif
		//xiaohui add 2015-03-24

		if(GlobalV::MY_RANK==0)
		{
//			upf.print_pseudo_upf( ofs );
			atoms[i].set_pseudo_nc( upf );

			log << "\n Read in pseudopotential file is " << pseudo_fn[i] << std::endl;
			ModuleBase::GlobalFunc::OUT(log,"pseudopotential type",atoms[i].pp_type);
			ModuleBase::GlobalFunc::OUT(log,"exchange-correlation functional", atoms[i].xc_func);
			ModuleBase::GlobalFunc::OUT(log,"nonlocal core correction", atoms[i].nlcc);
//			ModuleBase::GlobalFunc::OUT(log,"spin orbital",atoms[i].has_so);
			ModuleBase::GlobalFunc::OUT(log,"valence electrons", atoms[i].zv);
			ModuleBase::GlobalFunc::OUT(log,"lmax", atoms[i].lmax);
			ModuleBase::GlobalFunc::OUT(log,"number of zeta", atoms[i].nchi);
			ModuleBase::GlobalFunc::OUT(log,"number of projectors", atoms[i].nbeta);
			for(int ib=0; ib<atoms[i].nbeta; ib++)
			{
				ModuleBase::GlobalFunc::OUT(log,"L of projector", atoms[i].lll[ib]);
			}
//			ModuleBase::GlobalFunc::OUT(log,"Grid Mesh Number", atoms[i].mesh);
		}
		if(upf.functional_error == 1)
		{
			std::cout << "In Pseudopot_upf::read_pseudo_header : dft_functional from INPUT does not match that in pseudopot file" << std::endl;
			std::cout << "Please make sure this is what you need" << std::endl;
			atoms[i].xc_func = GlobalV::DFT_FUNCTIONAL;
			transform(atoms[i].xc_func.begin(), atoms[i].xc_func.end(), atoms[i].xc_func.begin(), (::toupper));
			if(GlobalV::MY_RANK==0)
			{
				log << "\n In Pseudopot_upf::read_pseudo_header : dft_functional from INPUT does not match that in pseudopot file" << std::endl;
				log << " Please make sure this is what you need" << std::endl;
				log << " XC functional updated to : " << GlobalV::DFT_FUNCTIONAL << std::endl;
				ModuleBase::GlobalFunc::OUT(log,"exchange-correlation functional", atoms[i].xc_func);
			}
		}
			
		//atoms[i].print_pseudo_us(ofs);
	}

//	if(GlobalV::MY_RANK==0)
//	{
//		ofs.close();
//	}
	return;
}


void UnitCell_pseudo::print_unitcell_pseudo(const std::string &fn)
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("UnitCell_pseudo","print_unitcell_pseudo");
	std::ofstream ofs( fn.c_str() );

	this->print_cell(ofs);
	for (int i = 0;i < ntype;i++)
	{
		atoms[i].print_Atom(ofs);
	}

	ofs.close();
	return;
}


#ifdef __MPI
void UnitCell_pseudo::bcast_unitcell_pseudo(void)
{
	Parallel_Common::bcast_int( meshx );
	Parallel_Common::bcast_int( natomwfc );
	Parallel_Common::bcast_int( lmax );
	Parallel_Common::bcast_int( lmax_ppwf );
	//Parallel_Common::bcast_double( GlobalC::CHR.nelec );

	bcast_unitcell();
}

void UnitCell_pseudo::bcast_unitcell_pseudo2(void)
{
	bcast_unitcell2();
}
#endif
