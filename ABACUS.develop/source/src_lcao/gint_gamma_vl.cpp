#include "gint_gamma.h"
#include "grid_technique.h"
#include "../module_ORB/ORB_read.h"
#include "../src_pw/global.h"
#include "src_global/blas_connector.h"

#include "global_fp.h" // mohan add 2021-01-30
#include "../src_global/ylm.h"
//#include <vector>

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
}

inline void cal_psir_ylm(
	const int na_grid,  // how many atoms on this (i,j,k) grid
	const int grid_index, // 1d index of FFT index (i,j,k)
	const double delta_r, // delta_r of the uniform FFT grid
	double phi, // radial wave functions
	const int*const colidx,  // count total number of atomis orbitals
	const int*const bsize,  // ??
	double*const*const psir_ylm, // bxyz * LD_pool
	int*const*const cal_flag) // whether the atom-grid distance is larger than cutoff
{
    for (int id=0; id<na_grid; id++)
    {
        // there are two parameters we want to know here:
        // in which bigcell of the meshball the atom is in?
        // what's the cartesian coordinate of the bigcell?
        const int mcell_index=GridT.bcell_start[grid_index] + id;

        const int iat=GridT.which_atom[mcell_index]; // index of atom
        const int it=ucell.iat2it[ iat ]; // index of atom type
        const Atom*const atom=&ucell.atoms[it];
		
        const int imcell=GridT.which_bigcell[mcell_index];

        // meshball_positions should be the bigcell position in meshball
        // to the center of meshball.
        // calculated in cartesian coordinates
        // the vector from the grid which is now being operated to the atom position.
        // in meshball language, is the vector from imcell to the center cel, plus
        // tau_in_bigcell.
		const double mt[3] = {
        	GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0],
        	GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1],
        	GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2]};

		// number of grids in each big cell (bxyz)
        for(int ib=0; ib<pw.bxyz; ib++)
        {
            double *p=&psir_ylm[ib][colidx[id]];
            // meshcell_pos: z is the fastest
			const double dr[3] = {
            	GridT.meshcell_pos[ib][0] + mt[0],
            	GridT.meshcell_pos[ib][1] + mt[1],
            	GridT.meshcell_pos[ib][2] + mt[2]};

            double distance = std::sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

            //if(distance[ib][id] > ORB.Phi[it].getRcut())
            if(distance > (ORB.Phi[it].getRcut()- 1.0e-15))
            {
				cal_flag[ib][id]=0;
                ZEROS(p, bsize[id]);
            }
			else
			{
				cal_flag[ib][id]=1;

				//if(distance[id] > GridT.orbital_rmax) continue;
				//    Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
				if (distance < 1.0E-9) distance += 1.0E-9;

				//------------------------------------------------------
				// spherical harmonic functions Ylm
				//------------------------------------------------------
				std::vector<double> ylma;
				Ylm::sph_harm (    ucell.atoms[it].nwl,
						dr[0] / distance,
						dr[1] / distance,
						dr[2] / distance,
						ylma);
				// these parameters are related to interpolation
				// because once the distance from atom to grid point is known,
				// we can obtain the parameters for interpolation and
				// store them first! these operations save lots of efforts.
				const double position=distance / delta_r;

				const int ip = static_cast<int>(position);
				const double dx = position - ip;
				const double dx2 = dx * dx;
				const double dx3 = dx2 * dx;

				const double c3 = 3.0*dx2-2.0*dx3;
				const double c1 = 1.0-c3;
				const double c2 = (dx-2.0*dx2+dx3)*delta_r;
				const double c4 = (dx3-dx2)*delta_r;

				for (int iw=0; iw< atom->nw; ++iw, ++p)
				{
					if ( atom->iw2_new[iw] )
					{
						const Numerical_Orbital_Lm &philn = ORB.Phi[it].PhiLN(
								atom->iw2l[iw],
								atom->iw2n[iw]);
						phi = c1*philn.psi_uniform[ip] + c2*philn.dpsi_uniform[ip]			 // radial wave functions
							+ c3*philn.psi_uniform[ip+1] + c4*philn.dpsi_uniform[ip+1];
					}
					*p=phi * ylma[atom->iw2_ylm[iw]];
				}
			}// end distance<=(ORB.Phi[it].getRcut()-1.0e-15)
        }// end ib
    }// end id
}

//inline void cal_meshball_vlocal(int na_grid, int LD_pool, int* block_iw, int* bsize, int* colidx,
void Gint_Gamma::cal_meshball_vlocal(
	const int na_grid,
	const int LD_pool,
	const int*const block_iw,
	const int*const bsize,
	const int*const colidx,
	const int*const*const cal_flag,
	const double*const vldr3,
	const double*const*const psir_ylm,
	double*const*const psir_vlbr3,
	const int lgd_now,
	double*const*const GridVlocal)
{
	const char transa='N', transb='T';
	const double alpha=1, beta=1;

	//int allnw=colidx[na_grid];
	for(int ib=0; ib<pw.bxyz; ++ib)
	{
        for(int ia=0; ia<na_grid; ++ia)
        {
            if(cal_flag[ib][ia]>0)
            {
                for(int i=colidx[ia]; i<colidx[ia+1]; ++i)
                {
                    psir_vlbr3[ib][i]=psir_ylm[ib][i]*vldr3[ib];
                }
            }
            else
            {
                for(int i=colidx[ia]; i<colidx[ia+1]; ++i)
                {
                    psir_vlbr3[ib][i]=0;
                }
            }

        }
	}

	for(int ia1=0; ia1<na_grid; ++ia1)
	{
		const int iw1_lo=block_iw[ia1];
		const int m=bsize[ia1];
		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int iw2_lo=block_iw[ia2];
			if(iw1_lo<=iw2_lo)
			{
                int first_ib=0, last_ib=0;
                for(int ib=0; ib<pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                    {
                        first_ib=ib;
                        break;
                    }
                }
                for(int ib=pw.bxyz-1; ib>=0; --ib)
                {
                    if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                    {
                        last_ib=ib+1;
                        break;
                    }
                }
                int ib_length=last_ib-first_ib;
                if(ib_length<=0) continue;

                int cal_pair_num=0;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    cal_pair_num+=cal_flag[ib][ia1]*cal_flag[ib][ia2];
                }

                int n=bsize[ia2];
//omp_set_lock(&lock);
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &n, &m, &ib_length, &alpha,
                        &psir_vlbr3[first_ib][colidx[ia2]], &LD_pool,
                        &psir_ylm[first_ib][colidx[ia1]], &LD_pool,
                        &beta, &GridVlocal[iw1_lo][iw2_lo], &lgd_now);
                }
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                        {
                            int k=1;
                            dgemm_(&transa, &transb, &n, &m, &k, &alpha,
                                &psir_vlbr3[ib][colidx[ia2]], &LD_pool,
                                &psir_ylm[ib][colidx[ia1]], &LD_pool,
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

    const int omp_threads = omp_get_max_threads();
	omp_set_num_threads(std::max(1,omp_threads/GridT.nbx));		// Peize Lin update 2021.01.20

#ifdef __OPENMP
	#pragma omp parallel
#endif
	{
		//OUT(ofs_running, "start calculate gamma_vlocal");

		// it's a uniform grid to save orbital values, so the delta_r is a constant.
		const double delta_r=ORB.dr_uniform;

		double phi=0.0;
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
			double *GridVlocal_pool=new double [lgd_now*lgd_now];
			ZEROS(GridVlocal_pool, lgd_now*lgd_now);

			double **GridVlocal_thread=new double*[lgd_now];
			for (int i=0; i<lgd_now; i++)
			{
				GridVlocal_thread[i]=&GridVlocal_pool[i*lgd_now];
			}
			Memory::record("Gint_Gamma","GridVlocal",lgd_now*lgd_now,"double");

			const int LD_pool = max_size*ucell.nwmax;

			double *psir_ylm_pool=new double[pw.bxyz*LD_pool];

			//------------------------------------------------------
			// atomic basis sets
			//------------------------------------------------------
			double **psir_ylm=new double *[pw.bxyz];
			for(int i=0; i<pw.bxyz; i++)
			{
				psir_ylm[i]=&psir_ylm_pool[i*LD_pool];
			}
			ZEROS(psir_ylm_pool, pw.bxyz*LD_pool);

			double *psir_vlbr3_pool=new double[pw.bxyz*LD_pool];
			double **psir_vlbr3=new double *[pw.bxyz];
			for(int i=0; i<pw.bxyz; i++)
			{
				psir_vlbr3[i]=&psir_vlbr3_pool[i*LD_pool];
			}
			ZEROS(psir_vlbr3_pool, pw.bxyz*LD_pool);

			//------------------------------------------------------
			// whether the atom-grid distance is larger than
			// cutoff
			//------------------------------------------------------
			int **cal_flag=new int*[pw.bxyz];
			for(int i=0; i<pw.bxyz; i++)
			{
				cal_flag[i]=new int[max_size];
			}

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
						double *vldr3 = get_vldr3(ncyz, ibx, jby, kbz);

						//------------------------------------------------------
						// index of wave functions for each block
						//------------------------------------------------------
						int *block_iw = get_block_iw(na_grid, grid_index, this->max_size);
						
						int* colidx = get_colidx(na_grid, grid_index);
						
						//------------------------------------------------------
						// band size: number of columns of a band
						//------------------------------------------------------------------
						int* bsize = get_bsize(na_grid, grid_index);
						
						//------------------------------------------------------------------
						// compute atomic basis phi(r) with both radial and angular parts
						//------------------------------------------------------------------
						cal_psir_ylm(na_grid, grid_index, delta_r, phi,
						colidx, bsize,  psir_ylm, cal_flag);

						//------------------------------------------------------------------
						// calculate <phi_i|V|phi_j>
						//------------------------------------------------------------------
						cal_meshball_vlocal(na_grid, LD_pool, block_iw, bsize, colidx, cal_flag,
						vldr3, psir_ylm, psir_vlbr3, lgd_now, GridVlocal_thread);
						
						free(vldr3);		vldr3=nullptr;
						free(block_iw);		block_iw=nullptr;
						free(colidx);		colidx=nullptr;
						free(bsize);		bsize=nullptr;
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
						GridVlocal[i][j] += GridVlocal_thread[i][j];
					}
				}
			}

			delete[] GridVlocal_thread;
			delete[] GridVlocal_pool;

			for(int i=0; i<pw.bxyz; i++)
			{
				delete[] cal_flag[i];
			}
			delete[] cal_flag;
			delete[] psir_vlbr3;
			delete[] psir_vlbr3_pool;
			delete[] psir_ylm;
			delete[] psir_ylm_pool;
		} // end of if(max_size>0 && lgd_now>0)
	} // end of #pragma omp parallel

    omp_set_num_threads(omp_threads);

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
