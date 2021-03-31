#include "local_orbital_wfc.h"
#include "../src_pw/global.h"
#include "../src_io/wf_local.h"

#include "global_fp.h" // mohan add 2021-01-30

Local_Orbital_wfc::Local_Orbital_wfc()
{
	allocate_flag = false;
	allocate_aug_flag = false;
	trace_aug = nullptr;	
	wfck_flag = false;	
	complex_flag = false;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{
	// used for force
	if(allocate_aug_flag)
	{
		if(GAMMA_ONLY_LOCAL)
		{
			for(int is=0; is<NSPIN; is++)
			{
				for(int i=0; i<NBANDS; i++) 
				{
					delete[] this->WFC_GAMMA_aug[is][i];
				}
				delete[] this->WFC_GAMMA_aug[is];
			}
			delete[] this->WFC_GAMMA_aug;
		}
		//mohan add 2012-01-09
		else
		{
			for(int ik=0; ik<kv.nks; ++ik)
			{
				for(int ib=0; ib<NBANDS; ++ib)
				{
					delete[] this->WFC_K_aug[ik][ib];
				}
				delete[] this->WFC_K_aug[ik];
			}
			delete[] this->WFC_K_aug;
		}
	}

	// used for k-points.
	if(complex_flag && this->wfck_flag)
	{
		for(int i=0; i<kv.nks; i++)
		{
			//for(int j=0; j<NBANDS; j++)
			//{
			//	delete[] this->WFC_K[i][j];
			//}
			delete[] this->WFC_K[i];
			//cout<<"delete WFC_K["<<i<<"] success"<<endl;
		}
		delete[] this->WFC_K;
		//cout<<"delete WFC_K success"<<endl;
		if(NLOCAL!= 0 )
		{
			delete[] this->WFC_K_POOL;
			//cout<<"delete WFC_K_POOL success"<<endl;
		}
	}

}

void Local_Orbital_wfc::allocate_k(const Grid_Technique &gt)
{
	TITLE("Local_Orbital_wfc","allocate_k");
	if(NLOCAL < NBANDS)
	{
		WARNING_QUIT("Local_Orbital_wfc::allocate","NLOCAL<NBANDS");
	}

	// mohan add the flag 2011-03-02
	// allocate the first part (only once!).
	if(this->wfck_flag == false)
	{
		this->WFC_K = new complex<double>**[kv.nks];
		for(int ik=0; ik<kv.nks; ik++)
		{
			this->WFC_K[ik] = new complex<double>*[NBANDS];
		}
		this->wfck_flag = true;
	}
	
	if(this->complex_flag)
	{
		delete[] this->WFC_K_POOL;
		this->complex_flag = false;
	}
	// allocate the second part.
	//if(gt.lgd != 0) xiaohui modify 2015-02-04, fixed memory bug
	//if(gt.lgd != 0 && this->complex_flag == false)
	if(gt.lgd != 0)
	{
		//cout<<"gt.lgd="<<gt.lgd<<" ; NLOCAL="<<NLOCAL<<endl; //delete 2015-09-06, xiaohui
		const int page=NBANDS*gt.lgd;
		this->WFC_K_POOL=new complex<double> [kv.nks*page];
		ZEROS(WFC_K_POOL, kv.nks*page);
		for(int ik=0; ik<kv.nks; ik++)
		{
			for(int ib=0; ib<NBANDS; ib++)
			{
				this->WFC_K[ik][ib] = &WFC_K_POOL[ik*page+ib*gt.lgd];
				//cout<<"ik="<<ik<<" ib="<<ib<<endl<<"WFC_K address: "<<WFC_K[ik][ib]<<" WFC_K_POOL address: "<<&WFC_K_POOL[ik*page+ib*gt.lgd]<<endl;
			}
			//cout<<"set WFC_K pointer success, ik: "<<ik<<endl;
			Memory::record("LocalOrbital_Coef","WFC_K",NBANDS*NLOCAL,"cdouble");
			//OUT(ofs_running,"MemoryForWaveFunctions (MB)",mem);
			//cout<<"WFC_K["<<ik<<"] use "<<mem<<" MB"<<endl;
			this->complex_flag = true;
		}
	}

	if(wf.start_wfc == "atomic" )
	{
		
	}
	else if(wf.start_wfc == "file")
	{
		int error;
		cout << " Read in wave functions files: " << kv.nkstot << endl;
		for(int ik=0; ik<kv.nkstot; ++ik)
		{
			ofs_running << " Read in wave functions " << ik + 1 << endl;
			error = WF_Local::read_lowf_complex( this->WFC_K[ik], ik );
		}
#ifdef __MPI
		Parallel_Common::bcast_int(error);
#endif
		ofs_running << " Error=" << error << endl;
		if(error==1)
		{
			WARNING_QUIT("Local_Orbital_wfc","Can't find the wave function file: LOWF.dat");
		}
		else if(error==2)
		{
			WARNING_QUIT("Local_Orbital_wfc","In wave function file, band number doesn't match");
		}
		else if(error==3)
		{
			WARNING_QUIT("Local_Orbital_wfc","In wave function file, nlocal doesn't match");
		}
		else if(error==4)
		{
			WARNING_QUIT("Local_Orbital_wfc","In k-dependent wave function file, k point is not correct");
		}
	}
	else
	{
		WARNING_QUIT("Local_Orbital_wfc","check the parameter: start_wfc");
	}

	return;
}


void Local_Orbital_wfc::set_trace_aug(const Grid_Technique &gt)
{
	TITLE("Local_Orbital_wfc","set_trace_aug");
	timer::tick("Local_Orbital_wfc","set_trace_aug",'D');
	// this function must be called after ParaO.trace_loc_row
	// , ParaO.trace_loc_col and GridT.trace_lo have been called.

	if(OUT_LEVEL != "m") 
	{
		ofs_running << "\n SETUP ARRAY FOR EXTRA WAVE FUNCTIONS" << endl;
	}

	bool* occ2d = new bool[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		occ2d[i] = false;
	}

	//------------------------------
	// 2d parallel of H, S matrix.
	// because in some terms of force,
	// we need to adapted to this. 
	//------------------------------
	for(int i=0; i<NLOCAL; i++)
	{
		const int mu = ParaO.trace_loc_row[i];
		const int nu = ParaO.trace_loc_col[i];
		if(mu>=0 || nu>=0)
		{
			occ2d[i] = true;
		}
	}

	//(1) init dimension of c_aug
	this->daug = 0;
	
	//(2) global positions of elementes in c_aug
	delete[] trace_aug;
	trace_aug = new int[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		// this -1 is important.
		trace_aug[i] = -1;
		if(occ2d[i])
		{
			if(gt.trace_lo[i]<0) 
			{
				this->trace_aug[i] = daug;
				++daug;
			}
		}
	}

	delete[] occ2d;

	//---------------------------------
	//second part: prepare for c_aug.
	//---------------------------------
	static bool first = true;
	if(first)
	{
		if(GAMMA_ONLY_LOCAL)
		{
			this->WFC_GAMMA_aug = new double**[NSPIN];
			for(int is=0; is<NSPIN; ++is)
			{
				this->WFC_GAMMA_aug[is] = new double*[NBANDS];
			}
		}
		else //mohan add 2012-01-08
		{
			this->WFC_K_aug = new complex<double>**[kv.nks];
			for(int ik=0; ik<kv.nks; ++ik)
			{
				this->WFC_K_aug[ik] = new complex<double>*[NBANDS];
			}
		}
		first=false;
	}	
	
	if(allocate_aug_flag)
	{
		if(GAMMA_ONLY_LOCAL)
		{
			for(int is=0; is<NSPIN; ++is)
			{
				for(int i=0; i<NBANDS; ++i)
				{
					delete[] WFC_GAMMA_aug[is][i];
				}
			}
		}
		else
		{
			for(int ik=0; ik<kv.nks; ++ik)
			{
				for(int i=0; i<NBANDS; ++i)
				{
					delete[] WFC_K_aug[ik][i];
				}
			}
		}
		allocate_aug_flag = false;
	}

	if(daug != 0)
	{
		//------------------------------------------------------
    	// Create wave functions(Coefficients) in local basis.
    	// Same as evc in plane wave basis.
		//------------------------------------------------------
		if(GAMMA_ONLY_LOCAL)
		{
			for(int is=0; is<NSPIN; ++is)
			{	
				for(int i=0; i<NBANDS; ++i)
				{
					this->WFC_GAMMA_aug[is][i] = new double[daug];
					ZEROS(this->WFC_GAMMA_aug[is][i], daug);
				}
			}
			Memory::record("LocalOrbital_Coef","WFC_GAMMA_aug",NSPIN*NBANDS*daug,"double");
		}
		else // mohan add 2012-01-08
		{
			for(int ik=0; ik<kv.nks; ++ik)
			{
				for(int i=0; i<NBANDS; ++i)
				{
					this->WFC_K_aug[ik][i] = new complex<double>[daug];
					ZEROS(this->WFC_K_aug[ik][i], daug);
				}
			}
		}
		allocate_aug_flag = true;
	}

	if(OUT_LEVEL != "m") 
	{
		OUT(ofs_running,"daug",daug);
	}

	timer::tick("Local_Orbital_wfc","set_trace_aug",'D');
	return; 
}
