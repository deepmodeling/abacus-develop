#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "src_global/blas_connector.h"

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
}

// setup buffer parameters for tranforming 2D block-cyclic distributed DM matrix 
inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
    //return (localIndex/nblk*nprocs+myproc)*nblk+localIndex%nblk;
}


inline int localIndex(int globalIndex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalIndex%(nblk*nprocs))/nblk);
    return int(globalIndex/(nblk*nprocs))*nblk+globalIndex%nblk;
}


int Local_Orbital_Charge::setAlltoallvParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk)
{
    OUT(ofs_running,"enter setAlltoallvParameter, nblk", nblk);
    timer::tick("LCAO_Charge","newDM_index",'F');
    // setup blacs parameters
    int nprows=0;	
	int npcols=0;
	int nprocs=0;
    int myprow=0;
	int mypcol=0;
	int myproc=0;

    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    Cblacs_pinfo(&myproc, &nprocs);
    // OUT(ofs_running,"nprocs",nprocs);


    // init data arrays
    delete[] sender_size_process;
    sender_size_process=new int[nprocs];
    delete[] sender_displacement_process;
    sender_displacement_process=new int[nprocs];

    // OUT(ofs_running,"lgd_now",lgd_now);
    
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

    // OUT(ofs_running,"nprows",nprows);
    // OUT(ofs_running,"npcols",npcols);

    for(int i=0; i<nprows; ++i)
    {
        nRow_in_proc[i]=0;
    }
    for(int i=0; i<npcols; ++i)
    {
        nCol_in_proc[i]=0;
    }

    // count the number of elements to be received from each process
    for(int iGlobal=0; iGlobal<NLOCAL; ++iGlobal)
    {
        int iLocalGrid=GridT.trace_lo[iGlobal];
        if(iLocalGrid>=0)
        {
            //trace_global[iLocalGrid]=iGlobal;
            int p;
            trace_2D_row[iLocalGrid]=localIndex(iGlobal, nblk, nprows, p);
            trace_2D_prow[iLocalGrid]=p;
            nRow_in_proc[trace_2D_prow[iLocalGrid]]++;
            trace_2D_col[iLocalGrid]=localIndex(iGlobal, nblk, npcols, p);
            trace_2D_pcol[iLocalGrid]=p;
            nCol_in_proc[trace_2D_pcol[iLocalGrid]]++;
        }
    }
    // OUT(ofs_running,"NLOCAL",NLOCAL);
    receiver_displacement_process[0]=0;
    // OUT(ofs_running,"receiver_displacement_process[0]",receiver_displacement_process[0]);
    for(int pnum=0; pnum<nprocs; ++pnum)
    {
        int prow, pcol;
        Cblacs_pcoord(blacs_ctxt, pnum, &prow, &pcol);
        receiver_size_process[pnum]=nRow_in_proc[prow]*nCol_in_proc[pcol];
        if(NEW_DM>1)
        {
            OUT(ofs_running,"pnum",pnum);
            OUT(ofs_running,"prow",prow);
            OUT(ofs_running,"pcol",pcol);
            OUT(ofs_running,"nRow_in_proc",nRow_in_proc[prow]);
            OUT(ofs_running,"nCol_in_proc",nCol_in_proc[pcol]);
        }
        if(pnum>0)
        {
            receiver_displacement_process[pnum]=receiver_displacement_process[pnum-1]+receiver_size_process[pnum-1];
        }
    }
    // OUT(ofs_running,"last receiver_size_process",receiver_size_process[nprocs-1]);
    
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
            int src_idx=src_row*NLOCAL+src_col; // leanding dimension is set to NLOCAL for all processes

            int src_pcol=trace_2D_pcol[j];
            int src_proc=Cblacs_pnum(blacs_ctxt, src_prow, src_pcol);

            receiver_2D_index[pos[src_proc]]=src_idx;
            receiver_local_index[pos[src_proc]]=i*lgd_now+j;
            ++pos[src_proc];
        }
    }
    // OUT(ofs_running,"last receiver_2D_index",receiver_2D_index[lgd_now*lgd_now-1]);
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
    
    // OUT(ofs_running,"last sender_size_process",sender_size_process[nprocs-1]);
    // setup sender buffer
    sender_size=sender_size_process[0];
    sender_displacement_process[0]=0;
    for(int i=1; i<nprocs; ++i)
    {
        sender_size+=sender_size_process[i];
        sender_displacement_process[i]=sender_displacement_process[i-1]+sender_size_process[i-1];
    }
    
    // OUT(ofs_running,"sender_size",sender_size);
    delete[] sender_2D_index;
    sender_2D_index=new int[sender_size];
    delete[] sender_buffer;
    sender_buffer=new double[sender_size];

    // send the index of the elements to be received via MPI_Alltoall
    MPI_Alltoallv(receiver_2D_index, receiver_size_process, receiver_displacement_process, MPI_INT,
                  sender_2D_index, sender_size_process, sender_displacement_process, MPI_INT, comm_2D);


    if(NEW_DM>1)
    {
        ofs_running<<"receiver_size is "<<receiver_size<<" ; receiver_size of each process is:\n";
        for(int i=0; i<nprocs; ++i)
        {
            ofs_running<<receiver_size_process[i]<<" ";
        }
        ofs_running<<endl;
        ofs_running<<"sender_size is "<<sender_size<<" ; sender_size of each process is:\n";
        for(int i=0; i<nprocs; ++i)
        {
            ofs_running<<sender_size_process[i]<<" ";
        }
        ofs_running<<endl;
    }
    // OUT(ofs_running,"last sender_2D_index",sender_2D_index[lgd_now*lgd_now-1]);
    delete[] receiver_2D_index;
    timer::tick("LCAO_Charge","newDM_index",'F');
    return 0;
}


// allocate density kernel may change once the ion
// positions change
void Local_Orbital_Charge::allocate_gamma(const Grid_Technique &gt)
{
     TITLE("Local_Orbital_Charge","allocate_gamma");

    // mohan fix serious bug 2010-09-06
    this->lgd_now = gt.lgd;
    //xiaohui add 'OUT_LEVEL' line, 2015-09-16
    if(OUT_LEVEL != "m") OUT(ofs_running,"lgd_last",lgd_last);
    if(OUT_LEVEL != "m") OUT(ofs_running,"lgd_now",lgd_now);

    // mohan add 2010-07-01
    if(this->init_DM)
    {
		assert(lgd_last > 0);
		for (int is=0; is<NSPIN; is++)
		{
			delete[] DM[is];
			delete[] DM_pool[is];
		}
		delete[] DM;
		delete[] DM_pool;
		init_DM = false;
    }

    assert(lgd_now <= NLOCAL);

    // mohan update 2010-09-06
    if(lgd_now > 0)
    {
		this->DM = new double**[NSPIN];
		this->DM_pool = new double *[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			this->DM_pool[is]=new double [lgd_now*lgd_now];
			ZEROS(DM_pool[is], lgd_now*lgd_now);
			this->DM[is] = new double*[lgd_now];

			for (int i=0; i<lgd_now; i++)
			{
				DM[is][i] = &DM_pool[is][i*lgd_now];
			}
			Memory::record("LocalOrbital_Charge","Density_Kernal",NSPIN*lgd_now*lgd_now,"double");
		}
		this->init_DM = true;
        this->lgd_last = lgd_now;
        //xiaohui add 'OUT_LEVEL', 2015-09-16
        if(OUT_LEVEL != "m") ofs_running << " allocate DM , the dimension is " << lgd_now << endl;
    }
    else if(lgd_now == 0)
    {
        this->init_DM = false;
    }
    else
    {
        WARNING_QUIT("Local_Orbital_Charge::allocate","lgd<0!Something Wrong!");
    }
    
    setAlltoallvParameter(ParaO.comm_2D, ParaO.blacs_ctxt, ParaO.nb);

	// Peize Lin test 2019-01-16
    wfc_dm_2d.init();

    return;
}


// calculate the grid distributed DM matrix from 2D block-cyclic distributed DM matrix
// transform dm_gamma[is].c to this->DM[is]
void Local_Orbital_Charge::cal_dk_gamma_from_2D(void)
{
    timer::tick("LCAO_Charge","dm_2dTOgrid",'F');
    OUT(ofs_running,"cal_dk_gamma_from_2D, NSPIN", NSPIN);

    for(int is=0; is<NSPIN; ++is)
    {
        if(NEW_DM>1)
        // outputDM( ParaO.blacs_ctxt, ParaO.nb);
        {
            // int myid;
            // MPI_Comm_rank(MPI_COMM_WORLD, &myid);
            // if(myid==0)
            // {
            //     ofs_running<<"DM[0][0:1][0:1] before send:"<<endl;
            //     ofs_running<<"DM(0,0)"<<wfc_dm_2d.dm_gamma[is](0,0)<<" ";
            //     ofs_running<<"DM(0,1)"<<wfc_dm_2d.dm_gamma[is](1,0)<<endl;
            //     ofs_running<<"DM(1,0)"<<wfc_dm_2d.dm_gamma[is](0,1)<<" ";
            //     ofs_running<<"DM(1,1)"<<wfc_dm_2d.dm_gamma[is](1,1)<<endl;
            // }
            ofs_running<<"2D block parameters:\n"<<"nblk: "<<ParaO.nb<<endl;
            ofs_running<<"DM in 2D format:\n_________________________________________\n";
            for(int i=0; i<wfc_dm_2d.dm_gamma[is].nr; ++i)
            {
                for(int j=0; j<wfc_dm_2d.dm_gamma[is].nc; ++j)
                {
                    ofs_running<<wfc_dm_2d.dm_gamma[is](i,j)<<" ";
                }
                ofs_running<<endl;
            }
            ofs_running<<"=========================================\n";
        }

        // put data from dm_gamma[is] to sender index
        int nNONZERO=0;
        for(int i=0; i<sender_size; ++i)
        {
            const int idx=sender_2D_index[i];
            const int icol=idx%NLOCAL;
            const int irow=(idx-icol)/NLOCAL;
            // sender_buffer[i]=wfc_dm_2d.dm_gamma[is](irow,icol);
            sender_buffer[i]=wfc_dm_2d.dm_gamma[is](icol,irow); // sender_buffer is clomun major, 
                                                                // so the row and column index should be switched
            if(sender_buffer[i]!=0) ++nNONZERO;
        }
        if(NEW_DM>1) 
        {
            OUT(ofs_running,"number of non-zero elements in sender_buffer",nNONZERO);
            OUT(ofs_running,"sender_size",sender_size);
            OUT(ofs_running,"last sender_buffer",sender_buffer[sender_size-1]);
        }
        // transform data via MPI_Alltoallv
        MPI_Alltoallv(sender_buffer, sender_size_process, sender_displacement_process, MPI_DOUBLE,
                      receiver_buffer, receiver_size_process, receiver_displacement_process, MPI_DOUBLE, ParaO.comm_2D);
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

        if(NEW_DM>1)
        {
            OUT(ofs_running,"number of non-zero elements in receiver_buffer",nNONZERO);
            OUT(ofs_running,"receiver_size",receiver_size);
            OUT(ofs_running,"last receiver_buffer",receiver_buffer[receiver_size-1]);
            // ofs_running<<"DM[0][0:1][0:1] after receiver:"<<endl;
            // int idx0=GridT.trace_lo[0];
            // int idx1=GridT.trace_lo[1];
            // if(idx0>=0)
            // {
            //     ofs_running<<"DM(0,0)"<<DM[0][idx0][idx0]<<" ";
            // }
            // if(idx0>=0 && idx1>=0)
            // {
            //     ofs_running<<"DM(0,1)"<<DM[0][idx0][idx1]<<endl;
            //     ofs_running<<"DM(1,0)"<<DM[0][idx1][idx0]<<" ";
            // }
            // if(idx1>=0)
            // {
            //     ofs_running<<"DM(1,1)"<<DM[0][idx1][idx1]<<endl;
            // }
            //ofs_running<<DM[0][0][0]<<" "<<DM[0][0][1]<<endl;
            //ofs_running<<DM[0][1][0]<<" "<<DM[0][1][1]<<endl;
            ofs_running<<"DM in local grid:\n_________________________________________\n";
            for(int i=0; i<NLOCAL; ++i)
            {
                int ii=GridT.trace_lo[i];
                if(ii < 0) continue;
                for(int j=0; j<NLOCAL; ++j)
                {
                    int jj=GridT.trace_lo[j];
                    if(jj<0) continue;
                    ofs_running<<DM[is][ii][jj]<<" ";
                }
                ofs_running<<endl;
            }
            ofs_running<<"=========================================\n";
        }
    }
    timer::tick("LCAO_Charge","dm_2dTOgrid",'F');
	return;
}

//--------------------------------------------------------------
void Local_Orbital_Charge::cal_dk_gamma(void)
{
    TITLE("Local_Orbital_Charge","cal_density_kernal");
    timer::tick("LocalOrbital_Charge","cal_dk_gamma",'F');

    assert(NSPIN==kv.nks);

#ifdef __MPI //2015-09-06, xiaohui
	#if EXX_DM==2
	if( Exx_Global::Hybrid_Type::HF==exx_lcao.info.hybrid_type 
		|| Exx_Global::Hybrid_Type::PBE0==exx_lcao.info.hybrid_type 
		|| Exx_Global::Hybrid_Type::HSE==exx_lcao.info.hybrid_type )
		exx_lcao.DM_para.clear_DMr();
	#endif

	// Peize Lin update 2018-07-02
	for(int is=0; is<NSPIN; ++is )
	{
		for (int i=0; i<lgd_now; i++)
		{
			ZEROS(this->DM[is][i], lgd_now);
		}
	}

	// initialize
	int nprocs=0;
	int myid=0;
	//MPI_Status status;
	MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
	MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);


	// DSIZE: number of processors in diag world
	vector<int> bands_local(DSIZE);
	for (int id=0; id<DSIZE; id++)
	{
		bands_local[id] = (id<NBANDS%DSIZE) ? NBANDS/DSIZE+1 : NBANDS/DSIZE;
	}
	const int band_local = bands_local[DRANK];

	int lastband_in_proc = 0;
	for (int id=0, count_bands=0; id<DSIZE; id++)
	{
		count_bands += bands_local[id];
		if (count_bands >= NBANDS)
		{
			lastband_in_proc = id;
			break;
		}
	}

	matrix wg_local(NSPIN,band_local);
	for(int id=0, Total_Bands=0; id <= lastband_in_proc; ++id)
	{
		if(myid == id)
		{
			for(int is=0; is<NSPIN; is++)
			{
				for (int ib=0; ib<bands_local[myid]; ib++)
				{
					wg_local(is,ib) = wf.wg(is,Total_Bands+ib);
				}
			}
		}
		Total_Bands += bands_local[id];
	}

	for( int is=0; is<NSPIN; ++is )
	{
		matrix Z_wg( NLOCAL, band_local );
		if(myid <= lastband_in_proc)
		{
			for(int iw=0; iw<NLOCAL; iw++)
			{
				for(int ib=0; ib<bands_local[myid]; ib++)
				{
					Z_wg(iw,ib) = ParaO.Z_LOC[is][iw*bands_local[myid]+ib] * wg_local(is,ib);
				}
			}
		}

		const int row_col = (NLOCAL%300) ? NLOCAL/300+1 : NLOCAL/300;

		matrix Z_row;
		matrix Z_col;
		matrix rho_row_col;

		for(int row_count=0; row_count<row_col; row_count++)
		{
			const int row_remain = ( (row_count+1)*300 <= NLOCAL )
				? 300
				: NLOCAL - row_count*300;

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
				const int col_remain = ( (col_count+1)*300 <= NLOCAL )
					? 300
					: NLOCAL - col_count*300;

				Z_col.create( col_remain, band_local, false );
				for(int i_col=0; i_col<col_remain; i_col++)
				{
					const int col_index = i_col +col_count*300;
					for(int ib=0; ib<band_local; ib++)
					{
						Z_col(i_col,ib) = ParaO.Z_LOC[is][col_index*band_local+ib] ;
					}
				}

				rho_row_col.create( row_remain, col_remain, false );

				//for(int i_row=0; i_row<row_remain; i_row++)
				//  for(int i_col=0; i_col<col_remain; i_col++)
				//      for(int ib=0; ib<band_local; ib++)
				//          rho_row_col(i_row,i_col) += Z_row(i_row,ib) * Z_col(i_col,ib);

				LapackConnector::gemm(
						'N', 'T', 
						row_remain, col_remain, band_local,
						1, Z_row.c, band_local, Z_col.c, band_local,
						0, rho_row_col.c, col_remain);
				MPI_Barrier(DIAG_HPSEPS_WORLD);
				Parallel_Reduce::reduce_double_all( rho_row_col.c, row_remain*col_remain);

				if(GAMMA_ONLY_LOCAL)
				{
					for(int i_row=0; i_row<row_remain; i_row++)
					{
						const int row_index = row_count*300 + i_row;
						const int row_mu = GridT.trace_lo[row_index];
						if(row_mu<0)    continue;
						for(int i_col=0; i_col<col_remain; i_col++)
						{
							const int col_index = col_count*300 + i_col;
							const int col_nu = GridT.trace_lo[col_index];
							if(col_nu<0)    continue;
							this->DM[is][row_mu][col_nu] = rho_row_col(i_row,i_col);
						}
					}
				}
						
				#if EXX_DM==2
				if( Exx_Global::Hybrid_Type::HF==exx_lcao.info.hybrid_type 
					|| Exx_Global::Hybrid_Type::PBE0==exx_lcao.info.hybrid_type 
					|| Exx_Global::Hybrid_Type::HSE==exx_lcao.info.hybrid_type )
				{
					exx_lcao.DM_para.set_DM_gamma( rho_row_col, is, {row_count*300,col_count*300} );
				}
				#endif				
			}  // end for col_count
		}  // end for row_count

		ofs_running<<"DM[0][0:1][0:1] in cal_dk_gamma:"<<endl;

		int idx0=GridT.trace_lo[0];
		int idx1=GridT.trace_lo[1];

		if(idx0>=0)
		{
			ofs_running<<"DM(0,0)"<<DM[is][idx0][idx0]<<"\t";
		}
		if(idx0>=0 && idx1>=0)
		{
			ofs_running<<"DM(0,1)"<<DM[is][idx0][idx1]<<endl;
			ofs_running<<"DM(1,0)"<<DM[is][idx1][idx0]<<"\t";
		}
		if(idx1>=0)
		{
			ofs_running<<"DM(1,1)"<<DM[is][idx1][idx1]<<endl;
		}
	}  // end for is    
#endif //2015-09-06, xiaohui

    timer::tick("LocalOrbital_Charge","cal_dk_gamma",'F');
    return;
}
