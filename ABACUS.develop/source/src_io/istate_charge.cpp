#include "istate_charge.h"
#include "../src_pw/global.h"
#include "../src_pw/tools.h"

IState_Charge::IState_Charge(){}

IState_Charge::~IState_Charge(){}


void IState_Charge::begin(void)
{
	TITLE("IState_Charge","begin");

	cout << " Perform |psi(i)|^2 for selected bands." << endl;

	if(!GAMMA_ONLY_LOCAL)
	{
		WARNING_QUIT("IState_Charge::begin","Only available for GAMMA_ONLY_LOCAL now.");
	}

	int mode = 0;
	if(NBANDS_ISTATE > 0) mode = 1;
	else mode = 2;

	int fermi_band=0;
	int bands_below=0;
	int bands_above=0;

	// (2) cicle:
	// (2.1) calculate the selected density matrix
	// from wave functions.
	// (2.2) carry out the grid integration to
	// get the charge density.
	this->bands_picked = new int[NBANDS];
	ZEROS(bands_picked, NBANDS);

	// (1) 
	// (1.1) allocate the space for LOWF.WFC_GAMMA

	// (1.2) read in LOWF_GAMMA.dat
	OUT(ofs_running,"LOWF.allocate_flag",LOWF.get_allocate_flag());	
	cout << " number of electrons = " << CHR.nelec << endl;

	// mohan update 2011-03-21
	// if ucell is odd, it's correct,
	// if ucell is even, it's also correct.
	// +1.0e-8 in case like (2.999999999+1)/2
	fermi_band = static_cast<int>( (CHR.nelec+1)/2 + 1.0e-8 ) ;
	cout << " number of occupied bands = " << fermi_band << endl;

	if(mode == 1)
	{
		bands_below = NBANDS_ISTATE;
		bands_above = NBANDS_ISTATE;

		cout << " plot band decomposed charge density below fermi surface with " 
			<< bands_below << " bands." << endl;

		cout << " plot band decomposed charge density above fermi surface with " 
			<< bands_above << " bands." << endl;

		for(int ib=0; ib<NBANDS; ib++)
		{
			if( ib >= fermi_band - bands_below ) 
			{
				if( ib < fermi_band + bands_above)
				{
					bands_picked[ib] = 1;
				}
			}
		}
	}
	else if(mode == 2)
	{
		bool stop = false;
		stringstream ss;
		ss << global_out_dir << "istate.info";
		cout << " Open the file : " << ss.str() << endl; 
		if(MY_RANK==0)
		{
			ifstream ifs(ss.str().c_str());
			if(!ifs)
			{
				stop = true;
			}
			else
			{
				int band_index;
				for(int ib=0; ib<NBANDS; ++ib)
				{
					READ_VALUE(ifs, bands_picked[ib]);
				}
			}
		}

#ifdef __MPI
		Parallel_Common::bcast_bool(stop);
		Parallel_Common::bcast_int(bands_picked, NBANDS);
#endif
		if(stop)
		{
			ofs_warning << " Can't find the file : " << ss.str() << endl;
			WARNING_QUIT("IState_Charge::begin","can't find the istate file.");
		}
	}

	for(int ib=0; ib<NBANDS; ib++)
	{
		if(bands_picked[ib])
		{
			cout << " Perform band decomposed charge density for band " << ib+1 << endl;
			// (1)
			// This has been done once in LOOP_ions.
			// but here we need to done for each band.
			LOC.allocate_gamma(GridT);	
			
			// (2) calculate the density matrix for a partuclar 
			// band, whenever it is occupied or not.
			this->idmatrix(ib);

			// (3) zero out of charge density array. 
			for(int is=0; is<NSPIN; is++)
			{
				ZEROS( CHR.rho[is], pw.nrxx );
			}
			
			// (4) calculate charge density for a particular 
			// band.
   			UHM.GG.cal_rho(LOC.DM);
			CHR.save_rho_before_sum_band(); //xiaohui add 2014-12-09
			stringstream ss;
			ss << global_out_dir << "BAND" << ib + 1 << "_CHG";
			// 0 means definitely output charge density.
			for(int is=0; is<NSPIN; is++)
			{
				bool for_plot = true;
				CHR.write_rho(CHR.rho_save[is], is, 0, ss.str(), 3, for_plot );
			}
		}
	}

	delete[] bands_picked;
	return;
}


void IState_Charge::idmatrix(const int &ib)
{
	TITLE("IState_Charge","idmatrix");
		
	for(int is=0; is<NSPIN; is++)
	{
		for (int i=0; i<NLOCAL; i++)
		{
			const int mu_local = GridT.trace_lo[i];
			if ( mu_local >= 0)
			{
				// set a pointer.
				double *alpha = LOC.DM[is][mu_local];
				for (int j=i; j<NLOCAL; j++)
				{
					const int nu_local = GridT.trace_lo[j];
					if ( nu_local >= 0)
					{
						//---------------------------------------------------------
						// LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
						// we need to fix this function in near future.
						// -- mohan add 2021-02-09
						//---------------------------------------------------------
						WARNING_QUIT("IState_Charge::idmatrix","need to update LOWF.WFC_GAMMA");
						// 2 stands for degeneracy.
						//alpha[nu_local] += 2.0 * LOWF.WFC_GAMMA[is][ib][mu_local] * LOWF.WFC_GAMMA[is][ib][nu_local];
					}
				}
			}
		}
	}
	return;
}

