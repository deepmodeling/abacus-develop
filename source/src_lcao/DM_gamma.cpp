#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../src_io/wf_local.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
}

#ifdef __MPI
int Local_Orbital_Charge::setAlltoallvParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk)
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"enter setAlltoallvParameter, nblk", nblk);
    ModuleBase::timer::tick("LCAO_Charge","newDM_index");
    // setup blacs parameters
    int nprows=0;	
	int npcols=0;
	int nprocs=0;
    int myprow=0;
	int mypcol=0;
	int myproc=0;

    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    Cblacs_pinfo(&myproc, &nprocs);
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nprocs",nprocs);


    // init data arrays
    delete[] sender_size_process;
    sender_size_process=new int[nprocs];
    delete[] sender_displacement_process;
    sender_displacement_process=new int[nprocs];
    GlobalV::ofs_running << "checkpoint 2" << std::endl;
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_now",lgd_now);
    
    receiver_size=lgd_now*lgd_now;
    receiver_size_process=new int[nprocs];
    delete[] receiver_displacement_process;
    receiver_displacement_process=new int[nprocs];
    delete[] receiver_local_index;
    receiver_local_index=new int[receiver_size];
    delete[] receiver_buffer;
    receiver_buffer=new double[receiver_size];
    
    int *trace_2D_row=new int[lgd_now];
    int *trace_2D_col=new int[lgd_now];
    int *trace_2D_prow=new int[lgd_now];
    int *trace_2D_pcol=new int[lgd_now];
    //int *trace_global=new int[lgd_now];

    int *nRow_in_proc=new int[nprows];
    int *nCol_in_proc=new int[npcols];

    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nprows",nprows);
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npcols",npcols);

    for(int i=0; i<nprows; ++i)
    {
        nRow_in_proc[i]=0;
    }
    for(int i=0; i<npcols; ++i)
    {
        nCol_in_proc[i]=0;
    }

    // count the number of elements to be received from each process
    for(int iGlobal=0; iGlobal<GlobalV::NLOCAL; ++iGlobal)
    {
        int iLocalGrid=GlobalC::GridT.trace_lo[iGlobal];
        if(iLocalGrid>=0)
        {
            //trace_global[iLocalGrid]=iGlobal;
            int p;
            trace_2D_row[iLocalGrid]=Local_Orbital_wfc::localIndex(iGlobal, nblk, nprows, p);
            trace_2D_prow[iLocalGrid]=p;
            nRow_in_proc[trace_2D_prow[iLocalGrid]]++;
            trace_2D_col[iLocalGrid]=Local_Orbital_wfc::localIndex(iGlobal, nblk, npcols, p);
            trace_2D_pcol[iLocalGrid]=p;
            nCol_in_proc[trace_2D_pcol[iLocalGrid]]++;
        }
    }
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NLOCAL",GlobalV::NLOCAL);
    receiver_displacement_process[0]=0;
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"receiver_displacement_process[0]",receiver_displacement_process[0]);
    for(int pnum=0; pnum<nprocs; ++pnum)
    {
        int prow, pcol;
        Cblacs_pcoord(blacs_ctxt, pnum, &prow, &pcol);
        receiver_size_process[pnum]=nRow_in_proc[prow]*nCol_in_proc[pcol];

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"pnum",pnum);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"prow",prow);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"pcol",pcol);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nRow_in_proc",nRow_in_proc[prow]);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nCol_in_proc",nCol_in_proc[pcol]);

        if(pnum>0)
        {
            receiver_displacement_process[pnum]=receiver_displacement_process[pnum-1]+receiver_size_process[pnum-1];
        }
    }
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last receiver_size_process",receiver_size_process[nprocs-1]);
    
    // build the index to be received
    int* pos=new int[nprocs];
    int *receiver_2D_index=new int[receiver_size];
    for(int i=0; i<nprocs; ++i)
    {
        pos[i]=receiver_displacement_process[i];
    }
    for(int i=0; i<lgd_now; ++i)
    {
        int src_row=trace_2D_row[i];
        int src_prow=trace_2D_prow[i];
        for(int j=0; j<lgd_now; ++j)
        {
            int src_col=trace_2D_col[j];
            int src_idx=src_row*GlobalV::NLOCAL+src_col; // leanding dimension is set to GlobalV::NLOCAL for all processes

            int src_pcol=trace_2D_pcol[j];
            int src_proc=Cblacs_pnum(blacs_ctxt, src_prow, src_pcol);

            receiver_2D_index[pos[src_proc]]=src_idx;
            receiver_local_index[pos[src_proc]]=i*lgd_now+j;
            ++pos[src_proc];
        }
    }
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last receiver_2D_index",receiver_2D_index[lgd_now*lgd_now-1]);
    delete[] pos;
    delete[] trace_2D_row;
    delete[] trace_2D_col;
    delete[] trace_2D_prow;
    delete[] trace_2D_pcol;
    //delete[] trace_global;
    delete[] nRow_in_proc;
    delete[] nCol_in_proc;
    
    // send number of elements to be sent via MPI_Alltoall
    MPI_Alltoall(receiver_size_process, 1, MPI_INT,
                 sender_size_process, 1, MPI_INT, comm_2D);
    
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last sender_size_process",sender_size_process[nprocs-1]);
    // setup sender buffer
    sender_size=sender_size_process[0];
    sender_displacement_process[0]=0;
    for(int i=1; i<nprocs; ++i)
    {
        sender_size+=sender_size_process[i];
        sender_displacement_process[i]=sender_displacement_process[i-1]+sender_size_process[i-1];
    }
    
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"sender_size",sender_size);
    delete[] sender_2D_index;
    sender_2D_index=new int[sender_size];
    delete[] sender_buffer;
    sender_buffer=new double[sender_size];

    // send the index of the elements to be received via MPI_Alltoall
    MPI_Alltoallv(receiver_2D_index, receiver_size_process, receiver_displacement_process, MPI_INT,
                  sender_2D_index, sender_size_process, sender_displacement_process, MPI_INT, comm_2D);


    GlobalV::ofs_running << "receiver_size is " << receiver_size << " ; receiver_size of each process is:\n";
    for(int i=0; i<nprocs; ++i)
    {
        GlobalV::ofs_running<<receiver_size_process[i]<<" ";
    }
    GlobalV::ofs_running<<std::endl;
    GlobalV::ofs_running<<"sender_size is "<<sender_size<<" ; sender_size of each process is:\n";
    for(int i=0; i<nprocs; ++i)
    {
        GlobalV::ofs_running<<sender_size_process[i]<<" ";
    }
    GlobalV::ofs_running << std::endl;
        
        // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last sender_2D_index",sender_2D_index[lgd_now*lgd_now-1]);
    delete[] receiver_2D_index;
    ModuleBase::timer::tick("LCAO_Charge","newDM_index");
    return 0;
}
#endif

// allocate density kernel may change once the ion
// positions change
void Local_Orbital_Charge::allocate_gamma(const Grid_Technique &gt)
{
     ModuleBase::TITLE("Local_Orbital_Charge","allocate_gamma");

    // mohan fix serious bug 2010-09-06
    this->lgd_now = gt.lgd;
    //xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
    if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_last",lgd_last);
    if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_now",lgd_now);

    // mohan add 2010-07-01
    if(this->init_DM)
    {
		assert(lgd_last > 0);
		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] DM[is];
			delete[] DM_pool[is];
		}
		delete[] DM;
		delete[] DM_pool;
		init_DM = false;
    }

    assert(lgd_now <= GlobalV::NLOCAL);

    // mohan update 2010-09-06
    if(lgd_now > 0)
    {
		this->DM = new double**[GlobalV::NSPIN];
		this->DM_pool = new double *[GlobalV::NSPIN];
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->DM_pool[is]=new double [lgd_now*lgd_now];
			ModuleBase::GlobalFunc::ZEROS(DM_pool[is], lgd_now*lgd_now);
			this->DM[is] = new double*[lgd_now];

			for (int i=0; i<lgd_now; i++)
			{
				DM[is][i] = &DM_pool[is][i*lgd_now];
			}
			ModuleBase::Memory::record("LocalOrbital_Charge","Density_Kernal",GlobalV::NSPIN*lgd_now*lgd_now,"double");
		}
		this->init_DM = true;
        this->lgd_last = lgd_now;
        //xiaohui add 'GlobalV::OUT_LEVEL', 2015-09-16
        if(GlobalV::OUT_LEVEL != "m") GlobalV::ofs_running << " allocate DM , the dimension is " << lgd_now << std::endl;
    }
    else if(lgd_now == 0)
    {
        this->init_DM = false;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Local_Orbital_Charge::allocate","lgd<0!Something Wrong!");
    }
    
#ifdef __MPI
    setAlltoallvParameter(this->ParaV->comm_2D, this->ParaV->blacs_ctxt, this->ParaV->nb);
#endif

	// Peize Lin test 2019-01-16
    this->init_dm_2d();

	if(GlobalC::wf.init_wfc=="file")
	{
		this->gamma_file(gt, *this->LOWF);
	}

    return;
}

void Local_Orbital_Charge::gamma_file(const Grid_Technique& gt,
    Local_Orbital_wfc &lowf)
{
	ModuleBase::TITLE("Local_Orbital_Charge","gamma_file");

	int error;
	std::cout << " Read in gamma point wave function files " << std::endl;

	double **ctot;

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{

		lowf.wfc_gamma[is].create(this->ParaV->ncol, this->ParaV->nrow);
		lowf.wfc_gamma[is].zero_out();

		GlobalV::ofs_running << " Read in wave functions " << is << std::endl;
		error = WF_Local::read_lowf( ctot , is, lowf);
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

	}//loop ispin
}

void Local_Orbital_Charge::cal_dk_gamma_from_2D_pub(void)
{
    ModuleBase::TITLE("Local_Orbital_Charge","cal_dk_gamma_from_2D_pub");

	cal_dk_gamma_from_2D();
}
// calculate the grid distributed DM matrix from 2D block-cyclic distributed DM matrix
// transform dm_gamma[is].c to this->DM[is]
void Local_Orbital_Charge::cal_dk_gamma_from_2D(void)
{
    ModuleBase::timer::tick("LCAO_Charge","dm_2dTOgrid");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"cal_dk_gamma_from_2D, NSPIN", GlobalV::NSPIN);

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        // int myid;
        // MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        // if(myid==0)
        // {
        //     GlobalV::ofs_running<<"DM[0][0:1][0:1] before send:"<<std::endl;
        //     GlobalV::ofs_running<<"DM(0,0)"<<wfc_dm_2d.dm_gamma[is](0,0)<<" ";
        //     GlobalV::ofs_running<<"DM(0,1)"<<wfc_dm_2d.dm_gamma[is](1,0)<<std::endl;
        //     GlobalV::ofs_running<<"DM(1,0)"<<wfc_dm_2d.dm_gamma[is](0,1)<<" ";
        //     GlobalV::ofs_running<<"DM(1,1)"<<wfc_dm_2d.dm_gamma[is](1,1)<<std::endl;
        // }
        GlobalV::ofs_running<<"2D block parameters:\n"<<"nblk: "<<this->ParaV->nb<<std::endl;
        GlobalV::ofs_running<<"DM in 2D format:\n_________________________________________\n";
        for(int i=0; i<this->dm_gamma[is].nr; ++i)
        {
            for(int j=0; j<this->dm_gamma[is].nc; ++j)
            {
                GlobalV::ofs_running<<this->dm_gamma[is](i,j)<<" ";
            }
            GlobalV::ofs_running<<std::endl;
        }
        GlobalV::ofs_running<<"=========================================\n";

        // put data from dm_gamma[is] to sender index
        int nNONZERO=0;
        for(int i=0; i<sender_size; ++i)
        {
            const int idx=sender_2D_index[i];
            const int icol=idx%GlobalV::NLOCAL;
            const int irow=(idx-icol)/GlobalV::NLOCAL;
            // sender_buffer[i]=wfc_dm_2d.dm_gamma[is](irow,icol);
            sender_buffer[i]=this->dm_gamma[is](icol,irow); // sender_buffer is clomun major, 
                                                                // so the row and column index should be switched
            if(sender_buffer[i]!=0) ++nNONZERO;
        }

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of non-zero elements in sender_buffer",nNONZERO);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"sender_size",sender_size);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last sender_buffer",sender_buffer[sender_size-1]);

        // transform data via MPI_Alltoallv
        #ifdef __MPI
        MPI_Alltoallv(sender_buffer, sender_size_process, sender_displacement_process, MPI_DOUBLE,
                      receiver_buffer, receiver_size_process, receiver_displacement_process, MPI_DOUBLE, this->ParaV->comm_2D);
        #endif
        // put data from receiver buffer to this->DM[is]
        nNONZERO=0;
        // init DM[is]
        /*for(int i=0; i<lgd_now; ++i)
        {
            for(int j=0; j<lgd_now; ++j)
            {
                DM[is][i][j]=0;
            }
        }*/
        for(int i=0; i<receiver_size; ++i)
        {
            const int idx=receiver_local_index[i];
            const int icol=idx%lgd_now;
            const int irow=(idx-icol)/lgd_now;
            DM[is][irow][icol]=receiver_buffer[i];
            //DM[is][icol][irow]=receiver_buffer[i];
            if(receiver_buffer[i]!=0) ++nNONZERO;
        }


        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of non-zero elements in receiver_buffer",nNONZERO);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"receiver_size",receiver_size);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last receiver_buffer",receiver_buffer[receiver_size-1]);
        // GlobalV::ofs_running<<"DM[0][0:1][0:1] after receiver:"<<std::endl;
        // int idx0=GlobalC::GridT.trace_lo[0];
        // int idx1=GlobalC::GridT.trace_lo[1];
        // if(idx0>=0)
        // {
        //     GlobalV::ofs_running<<"DM(0,0)"<<DM[0][idx0][idx0]<<" ";
        // }
        // if(idx0>=0 && idx1>=0)
        // {
        //     GlobalV::ofs_running<<"DM(0,1)"<<DM[0][idx0][idx1]<<std::endl;
        //     GlobalV::ofs_running<<"DM(1,0)"<<DM[0][idx1][idx0]<<" ";
        // }
        // if(idx1>=0)
        // {
        //     GlobalV::ofs_running<<"DM(1,1)"<<DM[0][idx1][idx1]<<std::endl;
        // }
        //GlobalV::ofs_running<<DM[0][0][0]<<" "<<DM[0][0][1]<<std::endl;
        //GlobalV::ofs_running<<DM[0][1][0]<<" "<<DM[0][1][1]<<std::endl;
        GlobalV::ofs_running<<"DM in local grid:\n_________________________________________\n";
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            int ii=GlobalC::GridT.trace_lo[i];
            if(ii < 0) continue;
            for(int j=0; j<GlobalV::NLOCAL; ++j)
            {
                int jj=GlobalC::GridT.trace_lo[j];
                if(jj<0) continue;
                GlobalV::ofs_running<<DM[is][ii][jj]<<" ";
            }
            GlobalV::ofs_running<<std::endl;
        }
        GlobalV::ofs_running<<"=========================================\n";

    }
    ModuleBase::timer::tick("LCAO_Charge","dm_2dTOgrid");
	return;
}

//--------------------------------------------------------------
void Local_Orbital_Charge::cal_dk_gamma(void)
{
    ModuleBase::TITLE("Local_Orbital_Charge","cal_density_kernal");
    ModuleBase::timer::tick("LocalOrbital_Charge","cal_dk_gamma");

    assert(GlobalV::NSPIN==GlobalC::kv.nks);

#ifdef __MPI //2015-09-06, xiaohui
	#if EXX_DM==2
	if( Exx_Global::Hybrid_Type::HF==GlobalC::exx_lcao.info.hybrid_type 
		|| Exx_Global::Hybrid_Type::PBE0==GlobalC::exx_lcao.info.hybrid_type 
		|| Exx_Global::Hybrid_Type::HSE==GlobalC::exx_lcao.info.hybrid_type )
		GlobalC::exx_lcao.DM_para.clear_DMr();
	#endif

	// Peize Lin update 2018-07-02
	for(int is=0; is<GlobalV::NSPIN; ++is )
	{
		for (int i=0; i<lgd_now; i++)
		{
			ModuleBase::GlobalFunc::ZEROS(this->DM[is][i], lgd_now);
		}
	}

	// initialize
	int nprocs=0;
	int myid=0;
	//MPI_Status status;
	MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
	MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);


	// GlobalV::DSIZE: number of processors in diag world
	std::vector<int> bands_local(GlobalV::DSIZE);
	for (int id=0; id<GlobalV::DSIZE; id++)
	{
		bands_local[id] = (id<GlobalV::NBANDS%GlobalV::DSIZE) ? GlobalV::NBANDS/GlobalV::DSIZE+1 : GlobalV::NBANDS/GlobalV::DSIZE;
	}
	const int band_local = bands_local[GlobalV::DRANK];

	int lastband_in_proc = 0;
	for (int id=0, count_bands=0; id<GlobalV::DSIZE; id++)
	{
		count_bands += bands_local[id];
		if (count_bands >= GlobalV::NBANDS)
		{
			lastband_in_proc = id;
			break;
		}
	}

	ModuleBase::matrix wg_local(GlobalV::NSPIN,band_local);
	for(int id=0, Total_Bands=0; id <= lastband_in_proc; ++id)
	{
		if(myid == id)
		{
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				for (int ib=0; ib<bands_local[myid]; ib++)
				{
					wg_local(is,ib) = GlobalC::wf.wg(is,Total_Bands+ib);
				}
			}
		}
		Total_Bands += bands_local[id];
	}

	for( int is=0; is<GlobalV::NSPIN; ++is )
	{
		ModuleBase::matrix Z_wg( GlobalV::NLOCAL, band_local );
		if(myid <= lastband_in_proc)
		{
			for(int iw=0; iw<GlobalV::NLOCAL; iw++)
			{
				for(int ib=0; ib<bands_local[myid]; ib++)
				{
					Z_wg(iw,ib) = this->ParaV->Z_LOC[is][iw*bands_local[myid]+ib] * wg_local(is,ib);
				}
			}
		}

		const int row_col = (GlobalV::NLOCAL%300) ? GlobalV::NLOCAL/300+1 : GlobalV::NLOCAL/300;

		ModuleBase::matrix Z_row;
		ModuleBase::matrix Z_col;
		ModuleBase::matrix rho_row_col;

		for(int row_count=0; row_count<row_col; row_count++)
		{
			const int row_remain = ( (row_count+1)*300 <= GlobalV::NLOCAL )
				? 300
				: GlobalV::NLOCAL - row_count*300;

			Z_row.create( row_remain, band_local, false );
			for(int i_row=0; i_row<row_remain; i_row++)
			{
				const int row_index = row_count*300 + i_row;
				for(int ib=0; ib<band_local; ib++)
				{
					Z_row(i_row,ib) = Z_wg(row_index,ib);
				}
			}

			for(int col_count=0; col_count<row_col; col_count++)
			{
				const int col_remain = ( (col_count+1)*300 <= GlobalV::NLOCAL )
					? 300
					: GlobalV::NLOCAL - col_count*300;

				Z_col.create( col_remain, band_local, false );
				for(int i_col=0; i_col<col_remain; i_col++)
				{
					const int col_index = i_col +col_count*300;
					for(int ib=0; ib<band_local; ib++)
					{
						Z_col(i_col,ib) = this->ParaV->Z_LOC[is][col_index*band_local+ib] ;
					}
				}

				rho_row_col.create( row_remain, col_remain, false );

				//for(int i_row=0; i_row<row_remain; i_row++)
				//  for(int i_col=0; i_col<col_remain; i_col++)
				//      for(int ib=0; ib<band_local; ib++)
				//          rho_row_col(i_row,i_col) += Z_row(i_row,ib) * Z_col(i_col,ib);

				BlasConnector::gemm(
						'N', 'T', 
						row_remain, col_remain, band_local,
						1, Z_row.c, band_local, Z_col.c, band_local,
						0, rho_row_col.c, col_remain);
				MPI_Barrier(DIAG_HPSEPS_WORLD);
				Parallel_Reduce::reduce_double_all( rho_row_col.c, row_remain*col_remain);

				if(GlobalV::GAMMA_ONLY_LOCAL)
				{
					for(int i_row=0; i_row<row_remain; i_row++)
					{
						const int row_index = row_count*300 + i_row;
						const int row_mu = GlobalC::GridT.trace_lo[row_index];
						if(row_mu<0)    continue;
						for(int i_col=0; i_col<col_remain; i_col++)
						{
							const int col_index = col_count*300 + i_col;
							const int col_nu = GlobalC::GridT.trace_lo[col_index];
							if(col_nu<0)    continue;
							this->DM[is][row_mu][col_nu] = rho_row_col(i_row,i_col);
						}
					}
				}
						
				#if EXX_DM==2
				if( Exx_Global::Hybrid_Type::HF==GlobalC::exx_lcao.info.hybrid_type 
					|| Exx_Global::Hybrid_Type::PBE0==GlobalC::exx_lcao.info.hybrid_type 
					|| Exx_Global::Hybrid_Type::HSE==GlobalC::exx_lcao.info.hybrid_type )
				{
					GlobalC::exx_lcao.DM_para.set_DM_gamma( rho_row_col, is, {row_count*300,col_count*300} );
				}
				#endif				
			}  // end for col_count
		}  // end for row_count

		GlobalV::ofs_running<<"DM[0][0:1][0:1] in cal_dk_gamma:"<<std::endl;

		int idx0=GlobalC::GridT.trace_lo[0];
		int idx1=GlobalC::GridT.trace_lo[1];

		if(idx0>=0)
		{
			GlobalV::ofs_running<<"DM(0,0)"<<DM[is][idx0][idx0]<<"\t";
		}
		if(idx0>=0 && idx1>=0)
		{
			GlobalV::ofs_running<<"DM(0,1)"<<DM[is][idx0][idx1]<<std::endl;
			GlobalV::ofs_running<<"DM(1,0)"<<DM[is][idx1][idx0]<<"\t";
		}
		if(idx1>=0)
		{
			GlobalV::ofs_running<<"DM(1,1)"<<DM[is][idx1][idx1]<<std::endl;
		}
	}  // end for is    
#endif //2015-09-06, xiaohui

    ModuleBase::timer::tick("LocalOrbital_Charge","cal_dk_gamma");
    return;
}
