#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "module_ORB/ORB_read.h"
#include "src_pw/global.h"
#include "src_global/blas_connector.h"
#include <mkl_service.h>

#include "global_fp.h" // mohan add 2021-01-30

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
}

// atomic basis sets
// psir_vlbr3[pw.bxyz][LD_pool]
Gint_Tools::Array_Pool<double> get_psir_vlbr3(
	const int na_grid,
	const int LD_pool,
	const int*const block_index,
	const bool*const*const cal_flag,
	const double*const vldr3,
	const double*const*const psir_ylm)
{
	Gint_Tools::Array_Pool<double> psir_vlbr3(pw.bxyz, LD_pool);
	for(int ib=0; ib<pw.bxyz; ++ib)
	{
        for(int ia=0; ia<na_grid; ++ia)
        {
            if(cal_flag[ib][ia])
            {
                for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
                {
                    psir_vlbr3.ptr_2D[ib][i]=psir_ylm[ib][i]*vldr3[ib];
                }
            }
            else
            {
                for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
                {
                    psir_vlbr3.ptr_2D[ib][i]=0;
                }
            }

        }
	}
    return psir_vlbr3;
}

//inline void cal_meshball_vlocal(int na_grid, int LD_pool, int* block_iw, int* block_size, int* block_index,
void Gint_Gamma::cal_meshball_vlocal(
	const int na_grid,
	const int LD_pool,
	const int*const block_iw,
	const int*const block_size,
	const int*const block_index,
	const bool*const*const cal_flag,
	const double*const vldr3,
	const double*const*const psir_ylm,
	const double*const*const psir_vlbr3,
	const int lgd_now,
	double*const*const GridVlocal)
{
	const char transa='N', transb='T';
	const double alpha=1, beta=1;

	for(int ia1=0; ia1<na_grid; ++ia1)
	{
		const int iw1_lo=block_iw[ia1];
		const int m=block_size[ia1];
		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int iw2_lo=block_iw[ia2];
			if(iw1_lo<=iw2_lo)
			{
                int first_ib=0;
                for(int ib=0; ib<pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        first_ib=ib;
                        break;
                    }
                }
                int last_ib=0;
                for(int ib=pw.bxyz-1; ib>=0; --ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        last_ib=ib+1;
                        break;
                    }
                }
                const int ib_length = last_ib-first_ib;
                if(ib_length<=0) continue;

                int cal_pair_num=0;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
                }

                const int n=block_size[ia2];
//omp_set_lock(&lock);
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &n, &m, &ib_length, &alpha,
                        &psir_vlbr3[first_ib][block_index[ia2]], &LD_pool,
                        &psir_ylm[first_ib][block_index[ia1]], &LD_pool,
                        &beta, &GridVlocal[iw1_lo][iw2_lo], &lgd_now);
                }
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        {
                            int k=1;
                            dgemm_(&transa, &transb, &n, &m, &k, &alpha,
                                &psir_vlbr3[ib][block_index[ia2]], &LD_pool,
                                &psir_ylm[ib][block_index[ia1]], &LD_pool,
                                &beta, &GridVlocal[iw1_lo][iw2_lo], &lgd_now);
                        }
                    }
                }
//omp_unset_lock(&lock);

			}
		}
	}
}

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


//------------------------------------------------------------------
// mohan add notes: 2021-03-11
// this subroutine is used to transform data from grid integrals
// to 2D-block distribution
// s stands for 'sender' and r stands for 'receiver'
//------------------------------------------------------------------
inline int setBufferParameter(
	MPI_Comm comm_2D,
	int blacs_ctxt,
	int nblk,
	int& s_index_siz,
	int*& s_local_index,
	int*& s_siz_pro,
	int*& s_dis_pro,
	int& s_siz,
	double*& s_buffer,
	int& r_index_siz,
	int*& r_global_index,
	int*& r_siz_pro,
	int*& r_dis_pro,
	int& r_siz,
	double*& r_buffer)
{
	//-----------------------------------------
    // setup blacs parameters
	//-----------------------------------------
    int nprows, npcols, nprocs;
    int myprow, mypcol, myproc;

    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

	//-----------------------------------------
	// set index of current proor: myproc
	// set number of total proors: nprocs
	//-----------------------------------------
    Cblacs_pinfo(&myproc, &nprocs);

    // initialize data arrays
    delete[] s_siz_pro;
    delete[] s_dis_pro;
    delete[] r_siz_pro;
    delete[] r_dis_pro;

    s_siz_pro=new int[nprocs];
    s_dis_pro=new int[nprocs];
    r_siz_pro=new int[nprocs];
    r_dis_pro=new int[nprocs];

	//---------------------------------------------------------------------
    // build the local index to be sent to other pro (s_local_index),
    // the global index to be received from other pro (r_global_index),
    // the send/receive siz/dis for data exchange by MPI_Alltoall
	//---------------------------------------------------------------------
    s_index_siz=GridT.lgd*GridT.lgd*2;

    delete[] s_local_index;
    s_local_index=new int[s_index_siz];

    int *s_global_index=new int[s_index_siz];

    int pos=0;
    s_siz_pro[0]=0;
    for(int iproc=0; iproc<nprocs; ++iproc)
    {
        s_dis_pro[iproc]=pos;

        int iprow=0;
		int ipcol=0;
        Cblacs_pcoord(blacs_ctxt, iproc, &iprow, &ipcol);

        // find out the global index and local index of elements
		// in each pro based on 2D block cyclic distribution
        for(int irow=0, grow=0; grow<NLOCAL; ++irow)
        {
            grow=globalIndex(irow, nblk, nprows, iprow);
            int lrow=GridT.trace_lo[grow];

            if(lrow < 0 || grow >= NLOCAL) continue;

            for(int icol=0, gcol=0; gcol<NLOCAL; ++icol)
            {
                gcol=globalIndex(icol,nblk, npcols, ipcol);
                int lcol=GridT.trace_lo[gcol];
                if(lcol < 0 || gcol >= NLOCAL) continue;
                // if(pos<0 || pos >= current_s_index_siz)
                // {
                //     OUT(ofs_running, "pos error, pos:", pos);
                //     OUT(ofs_running, "irow:", irow);
                //     OUT(ofs_running, "icol:", icol);
                //     OUT(ofs_running, "grow:", grow);
                //     OUT(ofs_running, "gcol:", gcol);
                //     OUT(ofs_running, "lrow:", grow);
                //     OUT(ofs_running, "lcol:", gcol);
                // }
                s_global_index[pos]=grow;
                s_global_index[pos+1]=gcol;
                s_local_index[pos]=lrow;
                s_local_index[pos+1]=lcol;
                pos+=2;
            }
        }
        s_siz_pro[iproc]=pos-s_dis_pro[iproc];
    }

    MPI_Alltoall(s_siz_pro, 1, MPI_INT,
                 r_siz_pro, 1, MPI_INT, comm_2D);

    r_index_siz=r_siz_pro[0];
    r_dis_pro[0]=0;
    for(int i=1; i<nprocs; ++i)
    {
        r_index_siz+=r_siz_pro[i];
        r_dis_pro[i]=r_dis_pro[i-1]+r_siz_pro[i-1];
    }

	delete[] r_global_index;
	r_global_index=new int[r_index_siz];

    // send the global index in sendBuffer to recvBuffer
    MPI_Alltoallv(s_global_index, s_siz_pro, s_dis_pro, MPI_INT,
                  r_global_index, r_siz_pro, r_dis_pro, MPI_INT, comm_2D);

    delete [] s_global_index;

    // the s_siz_pro, s_dis_pro, r_siz_pro,
    // and r_dis_pro will be used in transfer s_buffer, which
    // is half siz of s_global_index
    // we have to rebuild the siz and dis for each pro
    for (int iproc=0; iproc < nprocs; ++iproc)
    {
        s_siz_pro[iproc]=s_siz_pro[iproc]/2;
        s_dis_pro[iproc]=s_dis_pro[iproc]/2;
        r_siz_pro[iproc]=r_siz_pro[iproc]/2;
        r_dis_pro[iproc]=r_dis_pro[iproc]/2;
    }

    s_siz=s_index_siz/2;
	delete[] s_buffer;
	s_buffer=new double[s_siz];

    r_siz=r_index_siz/2;
	delete[] r_buffer;
	r_buffer=new double[r_siz];

    return 0;
}



void Gint_Gamma::cal_vlocal(
    const double* vlocal_in)
{
    omp_init_lock(&lock);

    TITLE("Gint_Gamma","cal_vlocal");
    timer::tick("Gint_Gamma","cal_vlocal",'J');

    this->job=cal_local;
    this->vlocal=vlocal_in;
    this->save_atoms_on_grid(GridT);

    this->gamma_vlocal();

    omp_destroy_lock(&lock);

    timer::tick("Gint_Gamma","cal_vlocal",'J');
    return;
}



void Gint_Gamma::gamma_vlocal(void)						// Peize Lin update OpenMP 2020.09.27
{
    TITLE("Gint_Gamma","gamma_vlocal");
    timer::tick("Gint_Gamma","gamma_vlocal",'K');

    double ** GridVlocal = new double*[GridT.lgd];
    for (int i=0; i<GridT.lgd; i++)
    {
        GridVlocal[i] = new double[GridT.lgd];
        ZEROS(GridVlocal[i], GridT.lgd);
    }
			Memory::record("Gint_Gamma","GridVlocal",GridT.lgd*GridT.lgd,"double");

    const int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(std::max(1,mkl_threads/GridT.nbx));		// Peize Lin update 2021.01.20

#ifdef __OPENMP
	#pragma omp parallel
#endif
	{
		//OUT(ofs_running, "start calculate gamma_vlocal");

		// it's a uniform grid to save orbital values, so the delta_r is a constant.
		const double delta_r=ORB.dr_uniform;

		const int nbx=GridT.nbx;
		const int nby=GridT.nby;
		const int nbz_start=GridT.nbzp_start;
		const int nbz=GridT.nbzp;

		const int ncyz=pw.ncy*pw.nczp;

		const int lgd_now=GridT.lgd;
		if(max_size>0 && lgd_now>0)
		{
			//------------------------------------------------------
			// <phi | V_local | phi>
			//------------------------------------------------------
			Gint_Tools::Array_Pool<double> GridVlocal_thread(lgd_now, lgd_now);
			ZEROS(GridVlocal_thread.ptr_1D, lgd_now*lgd_now);
			Memory::record("Gint_Gamma","GridVlocal_therad",lgd_now*lgd_now,"double");

			const int LD_pool = max_size*ucell.nwmax;

#ifdef __OPENMP
			#pragma omp for
#endif
			for (int i=0; i< nbx; i++)
			{
				const int ibx=i*pw.bx;
				for (int j=0; j<nby; j++)
				{
					const int jby=j*pw.by;
					for (int k=nbz_start; k<nbz_start+nbz; k++) // FFT grid
					{
						const int grid_index = i*nby*nbz + j*nbz + (k-nbz_start);

						//------------------------------------------------------------------
						// get the value: how many atoms are involved in this grid (big cell)
						//------------------------------------------------------------------
						const int na_grid=GridT.how_many_atoms[ grid_index ];
						if(na_grid==0) continue;

						//------------------------------------------------------------------
						// kbz can be obtained using a previously stored array
						//------------------------------------------------------------------
						const int kbz=k*pw.bz-pw.nczp_start;

						//------------------------------------------------------------------
						// extract the local potentials.
						//------------------------------------------------------------------
						double *vldr3 = this->get_vldr3(ncyz, ibx, jby, kbz);

						//------------------------------------------------------
						// index of wave functions for each block
						//------------------------------------------------------
						int *block_iw = Gint_Tools::get_block_iw(na_grid, grid_index, this->max_size);
						
						int* block_index = Gint_Tools::get_block_index(na_grid, grid_index);
						
						//------------------------------------------------------
						// band size: number of columns of a band
						//------------------------------------------------------
						int* block_size = Gint_Tools::get_block_size(na_grid, grid_index);

						//------------------------------------------------------
						// whether the atom-grid distance is larger than cutoff
						//------------------------------------------------------
						bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);
						
						//------------------------------------------------------------------
						// compute atomic basis phi(r) with both radial and angular parts
						//------------------------------------------------------------------
						const Gint_Tools::Array_Pool<double> psir_ylm = Gint_Tools::cal_psir_ylm(
							na_grid, LD_pool, grid_index, delta_r,
							block_index, block_size, cal_flag);

                        const Gint_Tools::Array_Pool<double> psir_vlbr3 = get_psir_vlbr3(
                                na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);

						//------------------------------------------------------------------
						// calculate <phi_i|V|phi_j>
						//------------------------------------------------------------------
						this->cal_meshball_vlocal(
							na_grid, LD_pool, block_iw, block_size, block_index, cal_flag,
							vldr3, psir_ylm.ptr_2D, psir_vlbr3.ptr_2D, lgd_now, GridVlocal_thread.ptr_2D);
						
						free(vldr3);		vldr3=nullptr;
						free(block_iw);		block_iw=nullptr;
						free(block_index);		block_index=nullptr;
						free(block_size);		block_size=nullptr;

						for(int ib=0; ib<pw.bxyz; ++ib)
							free(cal_flag[ib]);
						free(cal_flag);			cal_flag=nullptr;
					}// k
				}// j
			}// i

#ifdef __OPENMP
			#pragma omp critical(cal_vl)
#endif
			{
				for (int i=0; i<lgd_now; i++)
				{
					for (int j=0; j<lgd_now; j++)
					{
						GridVlocal[i][j] += GridVlocal_thread.ptr_2D[i][j];
					}
				}
			}
		} // end of if(max_size>0 && lgd_now>0)
	} // end of #pragma omp parallel

    mkl_set_num_threads(mkl_threads);

    OUT(ofs_running, "temp variables are deleted");
    timer::tick("Gint_Gamma","gamma_vlocal",'K');
    MPI_Barrier(MPI_COMM_WORLD);
    timer::tick("Gint_Gamma","distri_vl",'K');

    // setup send buffer and receive buffer size
    // OUT(ofs_running, "Start transforming vlocal from grid distribute to 2D block");
    if(CHR.get_new_e_iteration())
    {
        timer::tick("Gint_Gamma","distri_vl_index",'K');
        setBufferParameter(ParaO.comm_2D, ParaO.blacs_ctxt, ParaO.nb,
                           ParaO.sender_index_size, ParaO.sender_local_index,
                           ParaO.sender_size_process, ParaO.sender_displacement_process,
                           ParaO.sender_size, ParaO.sender_buffer,
                           ParaO.receiver_index_size, ParaO.receiver_global_index,
                           ParaO.receiver_size_process, ParaO.receiver_displacement_process,
                           ParaO.receiver_size, ParaO.receiver_buffer);
        OUT(ofs_running, "vlocal exchange index is built");
        OUT(ofs_running, "buffer size(M):", (ParaO.sender_size+ParaO.receiver_size)*sizeof(double)/1024/1024);
        OUT(ofs_running, "buffer index size(M):", (ParaO.sender_index_size+ParaO.receiver_index_size)*sizeof(int)/1024/1024);
        timer::tick("Gint_Gamma","distri_vl_index",'K');
    }

    timer::tick("Gint_Gamma","distri_vl_value",'K');

    // put data to send buffer
    for(int i=0; i<ParaO.sender_index_size; i+=2)
    {
        const int irow=ParaO.sender_local_index[i];
        const int icol=ParaO.sender_local_index[i+1];
        if(irow<=icol)
		{
            ParaO.sender_buffer[i/2]=GridVlocal[irow][icol];
		}
        else
		{
            ParaO.sender_buffer[i/2]=GridVlocal[icol][irow];
		}
    }
    OUT(ofs_running, "vlocal data are put in sender_buffer, size(M):", ParaO.sender_size*8/1024/1024);

    // use mpi_alltoall to get local data
    MPI_Alltoallv(ParaO.sender_buffer, ParaO.sender_size_process, ParaO.sender_displacement_process, MPI_DOUBLE,
                  ParaO.receiver_buffer, ParaO.receiver_size_process,
					ParaO.receiver_displacement_process, MPI_DOUBLE, ParaO.comm_2D);

    OUT(ofs_running, "vlocal data are exchanged, received size(M):", ParaO.receiver_size*8/1024/1024);

    // put local data to H matrix
    for(int i=0; i<ParaO.receiver_index_size; i+=2)
    {
        const int g_row=ParaO.receiver_global_index[i];
        const int g_col=ParaO.receiver_global_index[i+1];
        // if(g_col<0 || g_col>=NLOCAL||g_row<0 || g_row>=NLOCAL)
        // {
        //     OUT(ofs_running, "index error, i:", i);
        //     OUT(ofs_running, "index：", ParaO.receiver_global_index[i]);
        //     OUT(ofs_running, "g_col:", g_col);
        //     OUT(ofs_running, "g_col:", g_col);
        // }
        LM.set_HSgamma(g_row,g_col,ParaO.receiver_buffer[i/2],'L');
    }

    timer::tick("Gint_Gamma","distri_vl_value",'K');
    timer::tick("Gint_Gamma","distri_vl",'K');

	for (int i=0; i<GridT.lgd; i++)
	{
		delete[] GridVlocal[i];
	}
	delete[] GridVlocal;

	//OUT(ofs_running, "ALL GridVlocal was calculated");
    return;
}
