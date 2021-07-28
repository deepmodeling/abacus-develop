#include "unitcell_pseudo.h"
#ifdef __LCAO
//#include "../module_orbital/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif

#ifndef __CELL
#include "../src_pw/global.h"
#endif

#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

//==========================================================
// Read pseudopotential according to the dir
//==========================================================
void UnitCell_pseudo::read_cell_pseudopots(const string &pp_dir)
{
	TITLE("UnitCell_pseudo","read_cell_pseudopots");
	// setup reading log for pseudopot_upf
	stringstream ss;
	ss << GlobalV::global_out_dir << "atom_pseudo.log";
	
	// Read in the atomic pseudo potentials
	string pp_address;
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
				error_ap = upf.average_p(); //added by zhengdy 2020-10-20
			}
		}

#ifdef __MPI
		Parallel_Common::bcast_int(error);
		Parallel_Common::bcast_int(error_ap);
#endif

		if(error_ap) 
		{
			WARNING_QUIT("UnitCell_pseudo::read_pseudopot","error when average the pseudopotential.");
		}

		if(error==1)
		{
			cout << " Pseudopotential directory now is : " << pp_address << endl;
			GlobalV::ofs_warning << " Pseudopotential directory now is : " << pp_address << endl;
			WARNING_QUIT("read_pseudopot","Couldn't find pseudopotential file.");
		}
		else if(error==2)
		{
			WARNING_QUIT("read_pseudopot","Pseudopotential data do not match.");
		}
		else if(error==3)
		{
			WARNING_QUIT("read_pseudopot","Check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n unreasonable large (should be near 1.0), ABACUS would quit. \n The solution is to turn off the wave functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
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

			GlobalV::ofs_running << "\n Read in pseudopotential file is " << pseudo_fn[i] << endl;
			OUT(GlobalV::ofs_running,"pseudopotential type",atoms[i].pp_type);
			OUT(GlobalV::ofs_running,"functional Ex", atoms[i].dft[0]);
			OUT(GlobalV::ofs_running,"functional Ec", atoms[i].dft[1]);
			OUT(GlobalV::ofs_running,"functional GCEx", atoms[i].dft[2]);
			OUT(GlobalV::ofs_running,"functional GCEc", atoms[i].dft[3]);
			OUT(GlobalV::ofs_running,"nonlocal core correction", atoms[i].nlcc);
//			OUT(GlobalV::ofs_running,"spin orbital",atoms[i].has_so);
			OUT(GlobalV::ofs_running,"valence electrons", atoms[i].zv);
			OUT(GlobalV::ofs_running,"lmax", atoms[i].lmax);
			OUT(GlobalV::ofs_running,"number of zeta", atoms[i].nchi);
			OUT(GlobalV::ofs_running,"number of projectors", atoms[i].nbeta);
			for(int ib=0; ib<atoms[i].nbeta; ib++)
			{
				OUT(GlobalV::ofs_running,"L of projector", atoms[i].lll[ib]);
			}
//			OUT(GlobalV::ofs_running,"Grid Mesh Number", atoms[i].mesh);
		}
		if(upf.functional_error == 1)
		{
			cout << "In Pseudopot_upf::read_pseudo_header : input xc functional does not match that in pseudopot file" << endl;
			cout << "Please make sure this is what you need" << endl;
			atoms[i].dft[0] = GlobalV::DFT_FUNCTIONAL;
			transform(atoms[i].dft[0].begin(), atoms[i].dft[0].end(), atoms[i].dft[0].begin(), (::toupper));
			atoms[i].dft[1].clear();
			atoms[i].dft[2].clear();
			atoms[i].dft[3].clear();
			if(GlobalV::MY_RANK==0)
			{
				GlobalV::ofs_running << "\n XC functional updated to : " << endl;
				OUT(GlobalV::ofs_running,"functional Ex", atoms[i].dft[0]);
				OUT(GlobalV::ofs_running,"functional Ec", atoms[i].dft[1]);
				OUT(GlobalV::ofs_running,"functional GCEx", atoms[i].dft[2]);
				OUT(GlobalV::ofs_running,"functional GCEc", atoms[i].dft[3]);
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


void UnitCell_pseudo::print_unitcell_pseudo(const string &fn, output &outp)
{
	if(GlobalV::test_pseudo_cell) TITLE("UnitCell_pseudo","print_unitcell_pseudo");
	ofstream ofs( fn.c_str() );

	this->print_cell(ofs, outp);
	for (int i = 0;i < ntype;i++)
	{
		atoms[i].print_Atom(ofs, outp);
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
