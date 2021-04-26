#include "gint_gamma.h"
#include "grid_technique.h"
#include "ORB_read.h"
#include "../src_pw/global.h"
#include "src_global/blas_connector.h"
#include <mkl_service.h>

#include "global_fp.h" // mohan add 2021-01-30
#include "../src_global/ylm.h"

void Gint_Gamma::setVindex(const int ncyz, const int ibx, const int jby, const int kbz, int* vindex) const
{
    int bindex = 0;
    // z is the fastest, 
    for(int ii=0; ii<pw.bx; ii++)
    {
        const int ipart = (ibx + ii) * ncyz + kbz;
        for(int jj=0; jj<pw.by; jj++)
        {
            const int jpart = (jby + jj) * pw.nczp + ipart;
            for(int kk=0; kk<pw.bz; kk++)
            {
                vindex[bindex] = kk + jpart; 
                ++bindex;
            }
        }
    }
}

void Gint_Gamma::cal_psir_ylm_rho(
	int size, // number of atoms on this grid 
	int grid_index,  
	double delta_r, // delta_r of uniform grid
	double** distance, 
	double* ylma, 
	int* at, 
	int* block_index, 
	int* block_iw, 
	int* block_size, 
	int** cal_flag, 
	double** psir_ylm)
{
	const Numerical_Orbital_Lm *pointer;
	double mt[3];
	double dr[3];
	block_index[0]=0;
	for (int id=0; id<size; id++)
	{
		// there are two parameters we want to know here:
		// in which bigcell of the meshball the atom in?
		// what's the cartesian coordinate of the bigcell?
		const int mcell_index=GridT.bcell_start[grid_index] + id;
		const int imcell=GridT.which_bigcell[mcell_index];

		const int iat=GridT.which_atom[mcell_index];
		at[id]=iat;
		
		const int it=ucell.iat2it[iat];
		const int ia=ucell.iat2ia[iat];
		const int start=ucell.itiaiw2iwt(it, ia, 0);
		block_iw[id]=GridT.trace_lo[start];
		Atom* atom=&ucell.atoms[it];
		block_size[id]=atom->nw;
		block_index[id+1]=block_index[id]+atom->nw;
		// meshball_positions should be the bigcell position in meshball
		// to the center of meshball.
		// calculated in cartesian coordinates
		// the vector from the grid which is now being operated to the atom position.
		// in meshball language, is the vector from imcell to the center cel, plus
		// tau_in_bigcell.
		mt[0]=GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
		mt[1]=GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
		mt[2]=GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

		for(int ib=0; ib<pw.bxyz; ib++)
		{
			double *p=&psir_ylm[ib][block_index[id]];
			// meshcell_pos: z is the fastest
			dr[0]=GridT.meshcell_pos[ib][0] + mt[0]; 
			dr[1]=GridT.meshcell_pos[ib][1] + mt[1]; 
			dr[2]=GridT.meshcell_pos[ib][2] + mt[2]; 	

			distance[ib][id]=std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if(distance[ib][id] > (ORB.Phi[it].getRcut()- 1.0e-15)) 
			{
				cal_flag[ib][id]=0;
				ZEROS(p, block_size[id]);
				continue;
			}

			cal_flag[ib][id]=1;
			
			//if(distance[id] > GridT.orbital_rmax) continue;
			//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
			
			Ylm::sph_harm (	ucell.atoms[it].nwl,
					dr[0] / distance[ib][id],
					dr[1] / distance[ib][id],
					dr[2] / distance[ib][id],
					ylma);
			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position=distance[ib][id] / delta_r;
			const int ip=static_cast<int>(position);
			const double dx = position - ip;
			const double dx2 = dx * dx;
			const double dx3 = dx2 * dx;

			const double c3 = 3.0*dx2-2.0*dx3;
			const double c1 = 1.0-c3;
			const double c2 = (dx-2.0*dx2+dx3)*delta_r;
			const double c4 = (dx3-dx2)*delta_r;

			double phi=0;
			for (int iw=0; iw< atom->nw; ++iw, ++p)
			{
				if ( atom->iw2_new[iw] )
				{
					pointer=&ORB.Phi[it].PhiLN(
							atom->iw2l[iw],
							atom->iw2n[iw]);
					phi=c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
				}
				*p=phi * ylma[atom->iw2_ylm[iw]];
			} // end iw
		}// end ib
	}// end id
	return;
}

// can be done by GPU
void Gint_Gamma::cal_band_rho(
	int size, 
	int LD_pool, 
	int* block_iw, 
	int* bsize, 
	int* colidx,
	int** cal_flag, 
	double ** psir_ylm, 
	double **psir_DM, 
	double* psir_DM_pool, 
	int* vindex)
{

    //parameters for dsymm, dgemm and ddot
    char side='L', uplo='U';
    char transa='N', transb='N';
    double alpha_symm=1, alpha_gemm=2, beta=1;    
    int inc=1;

    for(int is=0; is<NSPIN; ++is)
    {
        ZEROS(psir_DM_pool, pw.bxyz*LD_pool);

        for (int ia1=0; ia1<size; ++ia1)
        {
            int iw1_lo=block_iw[ia1];

            //ia1==ia2, diagonal part
            // find the first ib and last ib for non-zeros cal_flag
            int first_ib=0, last_ib=0;
            for(int ib=0; ib<pw.bxyz; ++ib)
            {
                if(cal_flag[ib][ia1]>0)
                {
                    first_ib=ib;
                    break;
                }
            }
            for(int ib=pw.bxyz-1; ib>=0; --ib)
            {
                if(cal_flag[ib][ia1]>0)
                {
                    last_ib=ib+1;
                    break;
                }
            }
            int ib_length=last_ib-first_ib;
            if(ib_length<=0) continue;

            int cal_num=0;
            for(int ib=first_ib; ib<last_ib; ++ib)
            {
                cal_num+=cal_flag[ib][ia1];
            }
            // if enough cal_flag is nonzero
            if(cal_num>ib_length/4)
            {
                dsymm_(&side, &uplo, &bsize[ia1], &ib_length, 
                    &alpha_symm, &LOC.DM[is][iw1_lo][iw1_lo], &GridT.lgd, 
                    &psir_ylm[first_ib][colidx[ia1]], &LD_pool, 
                    &beta, &psir_DM[first_ib][colidx[ia1]], &LD_pool);
            }
            else
            {
                // int k=1;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    if(cal_flag[ib][ia1]>0)
                    {
                        dsymv_(&uplo, &bsize[ia1],
                            &alpha_symm, &LOC.DM[is][iw1_lo][iw1_lo], &GridT.lgd,
                            &psir_ylm[ib][colidx[ia1]], &inc,
                            &beta, &psir_DM[ib][colidx[ia1]], &inc);
                    }
                }
            }
            
            //OUT(ofs_running, "diagonal part of psir_DM done");
            for (int ia2=ia1+1; ia2<size; ++ia2)
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
                int iw2_lo=block_iw[ia2];
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &bsize[ia2], &ib_length, &bsize[ia1], 
                        &alpha_gemm, &LOC.DM[is][iw1_lo][iw2_lo], &GridT.lgd, 
                        &psir_ylm[first_ib][colidx[ia1]], &LD_pool, 
                        &beta, &psir_DM[first_ib][colidx[ia2]], &LD_pool);
                }
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1]>0 && cal_flag[ib][ia2]>0)
                        {
                            dgemv_(&transa, &bsize[ia2], &bsize[ia1], 
                                &alpha_gemm, &LOC.DM[is][iw1_lo][iw2_lo], &GridT.lgd,
                                &psir_ylm[ib][colidx[ia1]], &inc,
                                &beta, &psir_DM[ib][colidx[ia2]], &inc);
                        }
                    }
                }
                //OUT(ofs_running, "upper triangle part of psir_DM done, atom2", ia2);
            }// ia2
        } // ia1
    
        double *rhop = CHR.rho[is];
        for(int ib=0; ib<pw.bxyz; ++ib)
        {
            double r = ddot_(&colidx[size], psir_ylm[ib], &inc, psir_DM[ib], &inc);
            const int grid = vindex[ib];
            rhop[ grid ] += r;
        }
    } // end is
	return;
}

double Gint_Gamma::cal_rho(void)
{
    TITLE("Gint_Gamma","cal_rho");
    timer::tick("Gint_Gamma","cal_rho",'F');

    this->job = cal_charge;
    this->save_atoms_on_grid(GridT);

    double ne = 0.0;

	// I guess Peize add this, mohan 2021-01-31
	omp_init_lock(&lock);
    ne = this->gamma_charge();
	omp_destroy_lock(&lock);

    timer::tick("Gint_Gamma","cal_rho",'F');
    return ne;
}



double Gint_Gamma::gamma_charge(void)					// Peize Lin update OpenMP 2020.09.28
{
    TITLE("Gint_Gamma","gamma_charge");
    timer::tick("Gint_Gamma","gamma_charge",'I');    
    double sum = 0.0;//LiuXh 2016-01-10

	if(max_size)
    {
        const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(std::max(1,mkl_threads/GridT.nbx));			// Peize Lin update 2021.01.20
		
#ifdef __OPENMP
		#pragma omp parallel
#endif
		{
			// it's a uniform grid to save orbital values, so the delta_r is a constant.
			const double delta_r = ORB.dr_uniform;
		
			//if(omp_get_thread_num() == 0) ofs_running<<__FILE__<<__LINE__<<endl;

			// allocate 1
			int nnnmax=0;
			for(int T=0; T<ucell.ntype; T++)
			{
				nnnmax = max(nnnmax, nnn[T]);
			}

			// set up band matrix psir_ylm and psir_DM
			const int LD_pool=max_size*ucell.nwmax;
		
			double** distance = new double*[pw.bxyz];	// distance between atom and grid: [bxyz, maxsize]
			for(int i=0; i<pw.bxyz; ++i)
			{
				distance[i] = new double[max_size];
				ZEROS(distance[i], max_size);
			}
			int *block_size=new int[max_size];		//band size: number of columns of a band
			int *block_index=new int[max_size+1];
			int *at=new int[max_size];

			double *psir_ylm_pool=new double[pw.bxyz*LD_pool];
			ZEROS(psir_ylm_pool, pw.bxyz*LD_pool);
			
			double **psir_ylm=new double *[pw.bxyz];

			for(int i=0; i<pw.bxyz; ++i)
			{
				psir_ylm[i] = &psir_ylm_pool[i*LD_pool];
			}

			double *psir_DM_pool=new double[pw.bxyz*LD_pool];
			ZEROS(psir_DM_pool, pw.bxyz*LD_pool);

			double **psir_DM=new double *[pw.bxyz];

			for(int i=0; i<pw.bxyz; ++i)
			{
				psir_DM[i] = &psir_DM_pool[i*LD_pool];
			}

			int **cal_flag=new int*[pw.bxyz];

			for(int i=0; i<pw.bxyz; ++i)
			{
				cal_flag[i]=new int[max_size];
			}

			int *block_iw=new int[max_size]; // index of wave functions of each block;
		
			double* ylma = new double[nnnmax]; // Ylm for each atom: [bxyz, nnnmax]
			ZEROS(ylma, nnnmax);
		
			double v1 = 0.0;
			int* vindex=new int[pw.bxyz];
			ZEROS(vindex, pw.bxyz);
		
			const int nbx = GridT.nbx;
			const int nby = GridT.nby;
			const int nbz_start = GridT.nbzp_start;
			const int nbz = GridT.nbzp;
		
			const int ncyz = pw.ncy*pw.nczp; // mohan add 2012-03-25

#ifdef __OPENMP
			#pragma omp for
#endif
			for (int i=0; i<nbx; i++)
			{
				const int ibx = i*pw.bx;
				for (int j=0; j<nby; j++)
				{
					const int jby = j*pw.by;
					for (int k=nbz_start; k<nbz_start+nbz; k++)
					{
						const int kbz = k*pw.bz-pw.nczp_start;
		
						const int grid_index_thread = (k-nbz_start) + j * nbz + i * nby * nbz;
		
						// get the value: how many atoms has orbital value on this grid.
						//const int size = GridT.how_many_atoms[ this->grid_index ];
						const int size = GridT.how_many_atoms[ grid_index_thread ];
						if(size==0) continue;
						this->setVindex(ncyz, ibx, jby, kbz, vindex);
						
						this->cal_psir_ylm_rho(size, grid_index_thread, delta_r, distance, ylma,
								at, block_index, block_iw, block_size, 
								cal_flag, psir_ylm);
						

						cal_band_rho(size, LD_pool, block_iw, block_size, block_index,
								cal_flag, psir_ylm, psir_DM, psir_DM_pool, vindex);
					}// k
				}// j
			}// i
			
			delete[] vindex;
			delete[] ylma;
			
			delete[] block_iw;
			//OUT(ofs_running, "block_iw deleted");
			for(int i=0; i<pw.bxyz; i++)
			{
				delete[] cal_flag[i];
				delete[] distance[i];
			}
			delete[] cal_flag;
			delete[] distance;
			delete[] psir_DM;
			delete[] psir_DM_pool;
			delete[] psir_ylm;
			delete[] psir_ylm_pool;
			delete[] block_index;
			delete[] block_size;
			delete[] at;
		} // end of #pragma omp parallel
		
        for(int is=0; is<NSPIN; is++)
		{
            for (int ir=0; ir<pw.nrxx; ir++)
			{
                sum += CHR.rho[is][ir];
			}
		}
			
        mkl_set_num_threads(mkl_threads);
    } // end of if(max_size)
        
//ENDandRETURN:
    if(OUT_LEVEL != "m") OUT(ofs_running, "sum", sum);

    timer::tick("Gint_Gamma","reduce_charge",'J');
#ifdef __MPI
    Parallel_Reduce::reduce_double_pool( sum );
#endif
    timer::tick("Gint_Gamma","reduce_charge",'J');
    OUT(ofs_warning,"charge density sumed from grid", sum * ucell.omega/ pw.ncxyz);

    double ne = sum * ucell.omega / pw.ncxyz;
    //xiaohui add 'OUT_LEVEL', 2015-09-16
    if(OUT_LEVEL != "m") OUT(ofs_running, "ne", ne);
    timer::tick("Gint_Gamma","gamma_charge",'I');

    return ne;
}

