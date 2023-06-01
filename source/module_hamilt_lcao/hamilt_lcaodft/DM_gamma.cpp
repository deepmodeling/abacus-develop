#include "local_orbital_charge.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/blas_connector.h"
#include "module_io/read_wfc_nao.h"
#include "module_base/parallel_reduce.h"
#include "module_base/parallel_common.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

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
    ModuleBase::timer::tick("LOC","Alltoall");
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

    size_t memory_sender = 0;
    size_t memory_receiver = 0;
    size_t memory_other = 0;
    // init data arrays
    delete[] sender_size_process;
    sender_size_process=new int[nprocs];
    delete[] sender_displacement_process;
    sender_displacement_process=new int[nprocs];
    //GlobalV::ofs_running << "checkpoint 2" << std::endl;
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_now",lgd_now);
    memory_sender += nprocs * 2 * sizeof(int);
    
    receiver_size=lgd_now*lgd_now;
    delete[] receiver_size_process;
    receiver_size_process=new int[nprocs];
    delete[] receiver_displacement_process;
    receiver_displacement_process=new int[nprocs];
    delete[] receiver_local_index;
    receiver_local_index=new int[receiver_size];
    delete[] receiver_buffer;
    receiver_buffer=new double[receiver_size];

    memory_receiver += nprocs * 2 * sizeof(int) + receiver_size * (sizeof(int) + sizeof(double));
    
    int *trace_2D_row=new int[lgd_now];
    int *trace_2D_col=new int[lgd_now];
    int *trace_2D_prow=new int[lgd_now];
    int *trace_2D_pcol=new int[lgd_now];
    //int *trace_global=new int[lgd_now];

    int *nRow_in_proc=new int[nprows];
    int *nCol_in_proc=new int[npcols];

    memory_other += (lgd_now * 4 + nprows + npcols) * sizeof(int);

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
        int iLocalGrid=this->LOWF->gridt->trace_lo[iGlobal];
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
    memory_other += nprocs * sizeof(int);
    memory_receiver += receiver_size * sizeof(int);

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
    memory_sender += sender_size * (sizeof(int) * sizeof(double));

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
    ModuleBase::Memory::record("LOC::A2A_receiv", memory_receiver);
    ModuleBase::Memory::record("LOC::A2A_sender", memory_sender);
    ModuleBase::Memory::record("LOC::A2A_other", memory_other);
    ModuleBase::timer::tick("LOC","Alltoall");
    return 0;
}
#endif

// allocate density kernel may change once the ion
// positions change
void Local_Orbital_Charge::allocate_gamma(
                const int& lgd, 
                psi::Psi<double>* psid, 
                elecstate::ElecState* pelec,
                const int& nks)
{
     ModuleBase::TITLE("Local_Orbital_Charge","allocate_gamma");

    // mohan fix serious bug 2010-09-06
    this->lgd_now = lgd;
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
		}
		this->init_DM = true;
        this->lgd_last = lgd_now;
        ModuleBase::Memory::record("LOC::DM", sizeof(double) * GlobalV::NSPIN*lgd_now*lgd_now);
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
    this->init_dm_2d(nks);

    if (INPUT.init_wfc == "file")
    {
        this->gamma_file(psid, this->LOWF[0], pelec);
    }
    return;
}

void Local_Orbital_Charge::gamma_file(psi::Psi<double>* psid, Local_Orbital_wfc &lowf, elecstate::ElecState* pelec)
{
	ModuleBase::TITLE("Local_Orbital_Charge","gamma_file");

	int error;
	std::cout << " Read in gamma point wave function files " << std::endl;

	double **ctot;

    //allocate psi
    int ncol = this->ParaV->ncol_bands;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="lapack_gvx" || GlobalV::KS_SOLVER == "scalapack_gvx"
#ifdef __CUSOLVER_LCAO
    ||GlobalV::KS_SOLVER=="cusolver"
#endif
    )
    {
        ncol = this->ParaV->ncol;
    }
    if(psid == nullptr)
    {
        ModuleBase::WARNING_QUIT("gamma_file", "psid should be allocated first!");
    }
    else
    {
        psid->resize(GlobalV::NSPIN, ncol, this->ParaV->nrow);
    }
    ModuleBase::GlobalFunc::ZEROS( psid->get_pointer(), psid->size() );

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{

		GlobalV::ofs_running << " Read in wave functions " << is << std::endl;
		error = ModuleIO::read_wfc_nao( ctot , is, this->ParaV, psid, pelec);
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
#ifdef __DEBUG
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"cal_dk_gamma_from_2D, NSPIN", GlobalV::NSPIN);
#endif

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
        /*GlobalV::ofs_running<<"2D block parameters:\n"<<"nblk: "<<this->ParaV->nb<<std::endl;
        GlobalV::ofs_running<<"DM in 2D format:\n_________________________________________\n";
        for(int i=0; i<this->dm_gamma[is].nr; ++i)
        {
            for(int j=0; j<this->dm_gamma[is].nc; ++j)
            {
                GlobalV::ofs_running<<this->dm_gamma[is](i,j)<<" ";
            }
            GlobalV::ofs_running<<std::endl;
        }
        GlobalV::ofs_running<<"=========================================\n";*/

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

#ifdef __DEBUG
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of non-zero elements in sender_buffer",nNONZERO);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"sender_size",sender_size);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last sender_buffer",sender_buffer[sender_size-1]);
#endif

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

#ifdef __DEBUG
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of non-zero elements in receiver_buffer",nNONZERO);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"receiver_size",receiver_size);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last receiver_buffer",receiver_buffer[receiver_size-1]);
#endif
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
        /*GlobalV::ofs_running<<"DM in local grid:\n_________________________________________\n";
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
        GlobalV::ofs_running<<"=========================================\n";*/

    }
    ModuleBase::timer::tick("LCAO_Charge","dm_2dTOgrid");
	return;
}
