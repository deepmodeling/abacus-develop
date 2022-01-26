#include "local_orbital_wfc.h"
#include "../src_pw/global.h"
#include "../src_io/wf_local.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

Local_Orbital_wfc::Local_Orbital_wfc()
{
	allocate_flag = false;
	allocate_aug_flag = false;
	allocate_aug_backup_flag = false;
	trace_aug = nullptr;	
	wfck_flag = false;	
	wfck_backup_flag = false;
	complex_flag = false;
	complex_backup_flag = false;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{
	// used for force
	if(allocate_aug_flag)
	{
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		{
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				for(int ib=0; ib<GlobalV::NBANDS; ++ib)
				{
					delete[] this->WFC_K_aug[ik][ib];
				}
				delete[] this->WFC_K_aug[ik];
			}
			delete[] this->WFC_K_aug;
		}
	}
	if( INPUT.tddft==1 && allocate_aug_backup_flag)
	{
		cout<<"delete WFC_K_aug_backup"<<endl;
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		{
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				for(int ib=0; ib<GlobalV::NBANDS; ++ib)
				{
					delete[] this->WFC_K_aug_backup[ik][ib];
				}
				delete[] this->WFC_K_aug_backup[ik];
			}
			delete[] this->WFC_K_aug_backup;
		}
	}
	if (INPUT.tddft==1 && complex_backup_flag && this->wfck_backup_flag )
	{
		cout<<"delete WFC_K_backup and WFC_K_POOL_backup"<<endl;
		for(int i=0; i<GlobalC::kv.nks; i++)
		{
			delete[] this->WFC_K_backup[i];
		}
		delete[] this->WFC_K_backup;
		if(GlobalV::NLOCAL!= 0 )
		{
			delete[] this->WFC_K_POOL_backup;
		}
	}
	// used for k-points.
	if(complex_flag && this->wfck_flag)
	{
		for(int i=0; i<GlobalC::kv.nks; i++)
		{
			//for(int j=0; j<GlobalV::NBANDS; j++)
			//{
			//	delete[] this->WFC_K[i][j];
			//}
			delete[] this->WFC_K[i];
			//std::cout<<"delete WFC_K["<<i<<"] success"<<std::endl;
		}
		delete[] this->WFC_K;
		//std::cout<<"delete WFC_K success"<<std::endl;
		if(GlobalV::NLOCAL!= 0 )
		{
			delete[] this->WFC_K_POOL;
			//std::cout<<"delete WFC_K_POOL success"<<std::endl;
		}
	}

}

void Local_Orbital_wfc::allocate_k(const Grid_Technique &gt)
{
	ModuleBase::TITLE("Local_Orbital_wfc","allocate_k");
	if(GlobalV::NLOCAL < GlobalV::NBANDS)
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate","NLOCAL<NBANDS");
	}

	// mohan add the flag 2011-03-02
	// allocate the first part (only once!).
	if(this->wfck_flag == false)
	{
		this->WFC_K = new std::complex<double>**[GlobalC::kv.nks];
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			this->WFC_K[ik] = new std::complex<double>*[GlobalV::NBANDS];
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
		//std::cout<<"gt.lgd="<<gt.lgd<<" ; GlobalV::NLOCAL="<<GlobalV::NLOCAL<<std::endl; //delete 2015-09-06, xiaohui
		const int page=GlobalV::NBANDS*gt.lgd;
		this->WFC_K_POOL=new std::complex<double> [GlobalC::kv.nks*page];
		ModuleBase::GlobalFunc::ZEROS(WFC_K_POOL, GlobalC::kv.nks*page);
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ib++)
			{
				this->WFC_K[ik][ib] = &WFC_K_POOL[ik*page+ib*gt.lgd];
				//std::cout<<"ik="<<ik<<" ib="<<ib<<std::endl<<"WFC_K address: "<<WFC_K[ik][ib]<<" WFC_K_POOL address: "<<&WFC_K_POOL[ik*page+ib*gt.lgd]<<std::endl;
			}
			//std::cout<<"set WFC_K pointer success, ik: "<<ik<<std::endl;
			ModuleBase::Memory::record("LocalOrbital_Coef","WFC_K",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
			//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MemoryForWaveFunctions (MB)",mem);
			//std::cout<<"WFC_K["<<ik<<"] use "<<mem<<" MB"<<std::endl;
			this->complex_flag = true;
		}
	}

	if(GlobalC::wf.start_wfc == "atomic" )
	{
		
	}
	else if(GlobalC::wf.start_wfc == "file")
	{
		int error;
		std::cout << " Read in wave functions files: " << GlobalC::kv.nkstot << std::endl;
		for(int ik=0; ik<GlobalC::kv.nkstot; ++ik)
		{
			GlobalV::ofs_running << " Read in wave functions " << ik + 1 << std::endl;
			error = WF_Local::read_lowf_complex( this->WFC_K[ik], ik , 0);
		}
#ifdef __MPI
		Parallel_Common::bcast_int(error);
#endif
		GlobalV::ofs_running << " Error=" << error << std::endl;
		if(error==1)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","Can't find the wave function file: LOWF.dat");
		}
		else if(error==2)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, band number doesn't match");
		}
		else if(error==3)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, nlocal doesn't match");
		}
		else if(error==4)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In k-dependent wave function file, k point is not correct");
		}
	}
	else
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc","check the parameter: start_wfc");
	}

	return;
}

void Local_Orbital_wfc::allocate_k_backup(const Grid_Technique &gt)
{
	ModuleBase::TITLE("Local_Orbital_wfc","allocate_k_backup");
	if(GlobalV::NLOCAL < GlobalV::NBANDS)
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate","NLOCAL<NBANDS");
	}
	
	if (INPUT.tddft==1 && this->wfck_backup_flag == false )
	{
		cout<<"new WFC_K_backup"<<endl;
		this->WFC_K_backup = new std::complex<double>**[GlobalC::kv.nks];
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			this->WFC_K_backup[ik] = new std::complex<double>*[GlobalV::NBANDS];
		}
		this->wfck_backup_flag = true;
	}
	
	if( INPUT.tddft==1 && this->complex_backup_flag)
	{
		cout<<"(complex_flag) delete WFC_K_POOL_backup"<<endl;
		 delete[] this->WFC_K_POOL_backup;
		this->complex_backup_flag = false;
	}
	if( INPUT.tddft==1 && gt.lgd != 0 )
		{
			cout<<"new WFC_K_POOL_backup"<<endl;
			const int page=GlobalV::NBANDS*gt.lgd;
			this->WFC_K_POOL_backup=new std::complex<double> [GlobalC::kv.nks*page];
			ModuleBase::GlobalFunc::ZEROS(WFC_K_POOL_backup, GlobalC::kv.nks*page);
			for(int ik=0; ik<GlobalC::kv.nks; ik++)
			{
				for(int ib=0; ib<GlobalV::NBANDS; ib++)
				{
					this->WFC_K_backup[ik][ib] = &WFC_K_POOL_backup[ik*page+ib*gt.lgd];
				}
			}
			this->complex_backup_flag = true ;
		}
	return;
}

void Local_Orbital_wfc::set_trace_aug(const Grid_Technique &gt)
{
	ModuleBase::TITLE("Local_Orbital_wfc","set_trace_aug");
	ModuleBase::timer::tick("Local_Orbital_wfc","set_trace_aug");
	// this function must be called after GlobalC::ParaO.trace_loc_row
	// , GlobalC::ParaO.trace_loc_col and GlobalC::GridT.trace_lo have been called.

	if(GlobalV::OUT_LEVEL != "m") 
	{
		GlobalV::ofs_running << "\n SETUP ARRAY FOR EXTRA WAVE FUNCTIONS" << std::endl;
	}

	bool* occ2d = new bool[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		occ2d[i] = false;
	}

	//------------------------------
	// 2d parallel of H, S matrix.
	// because in some terms of force,
	// we need to adapted to this. 
	//------------------------------
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		const int mu = GlobalC::ParaO.trace_loc_row[i];
		const int nu = GlobalC::ParaO.trace_loc_col[i];
		if(mu>=0 || nu>=0)
		{
			occ2d[i] = true;
		}
	}

	//(1) init dimension of c_aug
	this->daug = 0;
	
	//(2) global positions of elementes in c_aug
	delete[] trace_aug;
	trace_aug = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
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
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		/*{
			this->WFC_GAMMA_aug = new double**[GlobalV::NSPIN];
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				this->WFC_GAMMA_aug[is] = new double*[GlobalV::NBANDS];
			}
		}
		else //mohan add 2012-01-08*/
		{
			this->WFC_K_aug = new std::complex<double>**[GlobalC::kv.nks];
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				this->WFC_K_aug[ik] = new std::complex<double>*[GlobalV::NBANDS];
			}
		}
		first=false;
	}	
	
	if(allocate_aug_flag)
	{
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		/*{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					delete[] WFC_GAMMA_aug[is][i];
				}
			}
		}
		else*/
		{
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				for(int i=0; i<GlobalV::NBANDS; ++i)
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
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		/*{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{	
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					this->WFC_GAMMA_aug[is][i] = new double[daug];
					ModuleBase::GlobalFunc::ZEROS(this->WFC_GAMMA_aug[is][i], daug);
				}
			}
			ModuleBase::Memory::record("LocalOrbital_Coef","WFC_GAMMA_aug",GlobalV::NSPIN*GlobalV::NBANDS*daug,"double");
		}
		else // mohan add 2012-01-08*/
		{
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					this->WFC_K_aug[ik][i] = new std::complex<double>[daug];
					ModuleBase::GlobalFunc::ZEROS(this->WFC_K_aug[ik][i], daug);
				}
			}
		}
		allocate_aug_flag = true;
	}

	if(GlobalV::OUT_LEVEL != "m") 
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"daug",daug);
	}

	ModuleBase::timer::tick("Local_Orbital_wfc","set_trace_aug");
	return; 
}

void Local_Orbital_wfc::set_trace_aug_backup(const Grid_Technique &gt)
{
	ModuleBase::TITLE("Local_Orbital_wfc","set_trace_aug");
	ModuleBase::timer::tick("Local_Orbital_wfc","set_trace_aug");
	// this function must be called after GlobalC::ParaO.trace_loc_row
	// , GlobalC::ParaO.trace_loc_col and GlobalC::GridT.trace_lo have been called.

	if(GlobalV::OUT_LEVEL != "m") 
	{
		GlobalV::ofs_running << "\n SETUP ARRAY FOR EXTRA WAVE FUNCTIONS" << std::endl;
	}

	bool* occ2d = new bool[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		occ2d[i] = false;
	}

	//------------------------------
	// 2d parallel of H, S matrix.
	// because in some terms of force,
	// we need to adapted to this. 
	//------------------------------
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		const int mu = GlobalC::ParaO.trace_loc_row[i];
		const int nu = GlobalC::ParaO.trace_loc_col[i];
		if(mu>=0 || nu>=0)
		{
			occ2d[i] = true;
		}
	}

	//(1) init dimension of c_aug
	this->daug = 0;
	
	//(2) global positions of elementes in c_aug
	delete[] trace_aug;
	trace_aug = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
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
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		/*{
			this->WFC_GAMMA_aug = new double**[GlobalV::NSPIN];
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				this->WFC_GAMMA_aug[is] = new double*[GlobalV::NBANDS];
			}
		}
		else //mohan add 2012-01-08*/
		{
			cout<<"new WFC_K_aug_backup the first and the second dimensions"<<endl;
			///*
			if (INPUT.tddft==1) 
			{
				this->WFC_K_aug_backup = new std::complex<double>**[GlobalC::kv.nks];
				for(int ik=0; ik<GlobalC::kv.nks; ++ik)
				{
					this->WFC_K_aug_backup[ik] = new std::complex<double>*[GlobalV::NBANDS];
				}
			}
			//*/
		}
		first=false;
	}	
	
	if(allocate_aug_backup_flag)
	{
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		/*{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					delete[] WFC_GAMMA_aug[is][i];
				}
			}
		}
		else*/
		{
			cout<<"delete WFC_K_aug_backup the third dimension "<<endl;
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					if (INPUT.tddft==1) delete[] WFC_K_aug_backup[ik][i];
				}
			}
		}
		allocate_aug_backup_flag = false;
	}

	if(daug != 0)
	{
		//------------------------------------------------------
    	// Create wave functions(Coefficients) in local basis.
    	// Same as evc in plane wave basis.
		//------------------------------------------------------
		if(!GlobalV::GAMMA_ONLY_LOCAL)
		/*{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{	
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					this->WFC_GAMMA_aug[is][i] = new double[daug];
					ModuleBase::GlobalFunc::ZEROS(this->WFC_GAMMA_aug[is][i], daug);
				}
			}
			ModuleBase::Memory::record("LocalOrbital_Coef","WFC_GAMMA_aug",GlobalV::NSPIN*GlobalV::NBANDS*daug,"double");
		}
		else // mohan add 2012-01-08*/
		{
			cout<<"new WFC_K_aug_backup the third dimension"<<endl;
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				for(int i=0; i<GlobalV::NBANDS; ++i)
				{
					///*
					if (INPUT.tddft==1) 
					{
						this->WFC_K_aug_backup[ik][i] = new std::complex<double>[daug];
						ModuleBase::GlobalFunc::ZEROS(this->WFC_K_aug_backup[ik][i], daug);
						//this->WFC_K_aug_backup[ik][i] = new std::complex<double>[GlobalV::NLOCAL];
						//ModuleBase::GlobalFunc::ZEROS(this->WFC_K_aug_backup[ik][i], GlobalV::NLOCAL);
					}
					//*/
				}
			}
		}
		allocate_aug_backup_flag = true;
	}

	if(GlobalV::OUT_LEVEL != "m") 
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"daug",daug);
	}

	ModuleBase::timer::tick("Local_Orbital_wfc","set_trace_aug");
	return; 
}