#include "subgrid_oper.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_pw/global.h"

SubGrid_oper::SubGrid_oper()
{
	trace_lo_tot = new int[1];	
	lgd=0;
	allocate_totwfc = false;
}
SubGrid_oper::~SubGrid_oper()
{
	delete[] trace_lo_tot;
}

//--------------------------------------------------
// because only the DIAG_WORLD processors
// have augmented wave functions,
// others will not have, we use this 
// augmented wave functions to calculate
// the force, not related to the grid integration
// Only DIAG_WORLD use density matrix (2D) only.
//--------------------------------------------------

void SubGrid_oper::cal_totwfc()
{
	ModuleBase::TITLE("SubGrid_oper","cal_totwfc");

	//-----------------------------------------
	// combine the wave functions index in the 
	// 'Grid group'
	//-----------------------------------------
	int *occupy = new int[GlobalV::NLOCAL];
	ModuleBase::GlobalFunc::ZEROS(occupy, GlobalV::NLOCAL);	
	for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
	{
		if(GlobalC::GridT.trace_lo[iw] >= 0)
		{
			occupy[iw] = 1;
		}
	}		

	// reduce occupy and get the full occupations.
#ifdef __MPI
	// mohan 2012-02-23
	Parallel_Reduce::reduce_int_grid(occupy, GlobalV::NLOCAL);
#endif

	/*
	for(int i=0; i<GlobalV::NLOCAL; ++i)
	{
		GlobalV::ofs_running << " i=" << i << " occupy=" << occupy[i] << std::endl;
	}
	*/

	delete[] trace_lo_tot;
	trace_lo_tot = new int[GlobalV::NLOCAL];
	ModuleBase::GlobalFunc::ZEROS(trace_lo_tot, GlobalV::NLOCAL);


	int count = 0;
	for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
	{
		// mohan fix bug 2012-02-23, shoud > 0,
		// not >= 0.
		if(occupy[iw]>0)
		{
			trace_lo_tot[iw]=count;
			++count;
		}
		else
		{
			trace_lo_tot[iw]=-1;
		}
	}
	//-----------------------------------------
	// calculate the total lgd.
	//-----------------------------------------
	this->lgd = count;
	
	//------------
	// for test
	//------------
	/*
	GlobalV::ofs_running << " trace_lo_tot" << std::endl;
	for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
	{
		GlobalV::ofs_running << " iw=" << iw << " trace_lo_tot=" << trace_lo_tot[iw] 
		<< " trace_lo=" << GlobalC::GridT.trace_lo[iw] 
		<< " occupy=" << occupy[iw] 
		<< std::endl;
	}
	*/

	if(GlobalV::GRANK==0)
	{
          	//xiaohui add 'GlobalV::OUT_LEVEL', 2015-09-16
		if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"GlobalC::SGO.lgd",lgd);
	
		if(lgd==0)
		{
			// need to allocate the pointer,
			// because it would be some array's
			// entrance parameters.
			// mohan update 2021-02-12
			this->totwfc = new double**[1];
			for(int is=0; is<1; ++is)
			{
				this->totwfc[is] = new double*[GlobalV::NBANDS];
			}
			this->allocate_totwfc=false;
			return;
		}
	
		assert(this->allocate_totwfc==false);

		static bool first_time=true;
		if(first_time)
		{
//			GlobalC::LOWF.init_Cij(1);
			first_time=false;
		}


		this->totwfc = new double**[1];
		for(int is=0; is<1; ++is)
		{
			this->totwfc[is] = new double*[GlobalV::NBANDS];
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				this->totwfc[is][ib] = new double[lgd];	
				ModuleBase::GlobalFunc::ZEROS(totwfc[is][ib], lgd);

				// mohan update 2012-02-10
				//if(GlobalV::DIAGO_TYPE!="cg") xiaohui modify 2013-09-02
				if(GlobalV::KS_SOLVER!="cg") //xiaohui add 2013-09-02
				{	
					for(int i=0; i<GlobalV::NLOCAL; ++i)
					{
						if(occupy[i])
						{
							if(ib==i)
							{
								totwfc[is][ib][trace_lo_tot[i]] = 1.0;
							}
						}
					}
				}
				else //use cg method
				{
					// bug: dimension of WFC_GAMMA is smaller than totwfc
					for(int i=0; i<lgd; ++i)
					{
		// mohan comment out 2021-02-09
		//				this->totwfc[is][ib][i] = GlobalC::LOWF.WFC_GAMMA[GlobalV::CURRENT_SPIN][ib][i]; //mohan update 2012-02-07	
					}
				}
			}	
		}

		allocate_totwfc = true;
	}

	delete[] occupy; //mohan fix bug 2012-03-25

	return;
}


void SubGrid_oper::dis_subwfc()
{
	ModuleBase::TITLE("SubGrid_oper","dis_subwfc");

#ifdef __MPI

//	std::cout << " distribute the wave functions " << std::endl;

	//------------------------------------------
	// bcast the eigenvalues
	//------------------------------------------
	for(int ik=0; ik<GlobalC::kv.nks; ++ik)
	{
		MPI_Bcast(GlobalC::wf.ekb[ik], GlobalV::NBANDS, MPI_DOUBLE, 0, GRID_WORLD);
	}	

	MPI_Status status;

	for(int i=0; i<GlobalV::GSIZE; ++i)
	{
		if(GlobalV::GRANK==0)
		{
			if(i==0)
			{
				//---------------------------------------------
				// Transfer the data from totwfc to WFC_GAMMA.
				//---------------------------------------------
				for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
				{
					const int mu1 = GlobalC::GridT.trace_lo[iw];
					if(mu1 >= 0)
					{
						const int mu2 = this->trace_lo_tot[iw];

						for(int ib=0; ib<GlobalV::NBANDS; ++ib)
						{
		// mohan comment out 2021-02-09
		//					GlobalC::LOWF.WFC_GAMMA[GlobalV::CURRENT_SPIN][ib][mu1] = this->totwfc[0][ib][mu2];
						}//ib
					}//mu1>=0
				}//iw
			}//i
			else
			{
				int tag;
				// receive trace_lo2
				tag = i * 10;
				int* trace_lo2 = new int[GlobalV::NLOCAL];
				MPI_Recv(trace_lo2, GlobalV::NLOCAL, MPI_INT, i, tag, GRID_WORLD, &status);

/*
				GlobalV::ofs_running << " Proc " << i << std::endl;
				for(int i=0; i<GlobalV::NLOCAL; ++i)
				{
					GlobalV::ofs_running << std::setw(5) << i << std::setw(10) << trace_lo2[i] << std::endl;
				}
				*/

				// receive lgd2
				int lgd2 = 0;
				tag = i * 10 + 1;
				MPI_Recv(&lgd2, 1, MPI_INT, i, tag, GRID_WORLD, &status);

//				GlobalV::ofs_running << " receive=" << tag << std::endl;
				
				// send csend
				double* csend = new double[GlobalV::NBANDS*lgd2];
				ModuleBase::GlobalFunc::ZEROS(csend, GlobalV::NBANDS*lgd2);

				for (int iw=0; iw<GlobalV::NLOCAL; iw++)
				{
					const int mu1 = trace_lo2[iw];
					if (mu1>=0)
					{
						const int mu2 = this->trace_lo_tot[iw]; 
						if(mu2<0)
						{
							GlobalV::ofs_running << " iw = " << iw << std::endl;
							GlobalV::ofs_running << " GlobalC::GridT.trace_lo=" << mu1 << std::endl;
							GlobalV::ofs_running << " trace_lo_tot=" << mu2 << std::endl;
							assert(mu2>=0);
						}

						for (int ib=0; ib<GlobalV::NBANDS; ib++)
						{
							csend[mu1*GlobalV::NBANDS+ib] = this->totwfc[0][ib][mu2];
						}
					}
				}

				tag = i * 10 + 2;
//				GlobalV::ofs_running << " send=" << tag << std::endl;
				MPI_Send(csend,GlobalV::NBANDS*lgd2,MPI_DOUBLE,i,tag,GRID_WORLD);
//				GlobalV::ofs_running << " send done." << std::endl;

				delete[] csend;
				delete[] trace_lo2;
			}
		}//GlobalV::GRANK=0
		else if(i==GlobalV::GRANK)
		{
			int tag;
			// send trace_lo
			tag = GlobalV::GRANK * 10;
			MPI_Send(GlobalC::GridT.trace_lo, GlobalV::NLOCAL, MPI_INT, 0, tag, GRID_WORLD);

			//GlobalV::ofs_running << " send1." << std::endl;

			// send GlobalC::GridT.lgd
			tag = GlobalV::GRANK * 10 + 1;
			MPI_Send(&GlobalC::GridT.lgd, 1, MPI_INT, 0, tag, GRID_WORLD);

			//GlobalV::ofs_running << " send2." << std::endl;

			// receive c
			double* crecv = new double[GlobalV::NBANDS*GlobalC::GridT.lgd];
			ModuleBase::GlobalFunc::ZEROS(crecv, GlobalV::NBANDS*GlobalC::GridT.lgd);

			tag = GlobalV::GRANK * 10 + 2;
//			GlobalV::ofs_running << " receive=" << tag << std::endl;
			MPI_Recv(crecv, GlobalV::NBANDS*GlobalC::GridT.lgd, MPI_DOUBLE, 0, tag, GRID_WORLD, &status);
//			GlobalV::ofs_running << " receive done." << std::endl;


			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				for (int mu=0; mu<GlobalC::GridT.lgd; ++mu)
				{
		// mohan comment out 2021-02-09
		//			GlobalC::LOWF.WFC_GAMMA[GlobalV::CURRENT_SPIN][ib][mu] = crecv[mu*GlobalV::NBANDS+ib];
				}
			}

			delete[] crecv;
		}
		MPI_Barrier(GRID_WORLD);
	}//end i


	//-------------------
	// Test
	//-------------------
	/*
	GlobalV::ofs_running << " WFC " << " CURRENT_SPIN=" << GlobalV::CURRENT_SPIN << std::endl;
	for(int i=0; i<GlobalV::NBANDS; ++i)
	{
		for(int j=0; j<GlobalV::NLOCAL; ++j)
		{
			const int mu = GlobalC::GridT.trace_lo[j];
			if(mu>=0)
			{
				//				if( abs(GlobalC::LOWF.WFC_GAMMA[0][i][mu] > 1.0e-8) )
				GlobalV::ofs_running << std::setw(5) << i+1 << std::setw(8) << j+1 
					<< std::setw(15) << GlobalC::LOWF.WFC_GAMMA[GlobalV::CURRENT_SPIN][i][mu] << std::endl; 
			}
		}
	}
	*/

#else

// this is for serial version

	//---------------------------------------------
	// Transfer the data from totwfc to WFC_GAMMA.
	//---------------------------------------------
	for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
	{
		const int mu1 = GlobalC::GridT.trace_lo[iw];
		if(mu1 >= 0)
		{
			const int mu2 = this->trace_lo_tot[iw];

			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
		// mohan comment out 2021-02-09
//				GlobalC::LOWF.WFC_GAMMA[GlobalV::CURRENT_SPIN][ib][mu1] = this->totwfc[0][ib][mu2];
			}//ib
		}//mu1>=0
	}//iw


#endif //mohan fix bug 2012-02-04

	if(GlobalV::GRANK==0)
	{
		// 1 spin because it can not be calculated two spin
		// at the same time.
		for(int is=0; is<1; ++is)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				if(allocate_totwfc)
				{
					delete[] this->totwfc[is][ib];
				}
			}
			delete[] this->totwfc[is];
		}
		delete[] this->totwfc;
		this->allocate_totwfc = false;
	}
	return;
}

