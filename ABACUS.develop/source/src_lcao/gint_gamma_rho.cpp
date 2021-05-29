#include "gint_gamma.h"
#include "grid_technique.h"
#include "../module_ORB/ORB_read.h"
#include "../src_pw/global.h"
#include "src_global/blas_connector.h"
#include <mkl_service.h>

#include "global_fp.h" // mohan add 2021-01-30
#include "../src_global/ylm.h"

void Gint_Gamma::cal_psir_ylm_rho(
	const int na_grid, // number of atoms on this grid 
	const int grid_index,  
	const double delta_r, // delta_r of uniform grid
	const int*const block_index,
	const int*const block_size, 
	bool*const*const cal_flag, 
	double*const*const psir_ylm)
{
	for (int id=0; id<na_grid; id++)
	{
		// there are two parameters we want to know here:
		// in which bigcell of the meshball the atom in?
		// what's the cartesian coordinate of the bigcell?
		const int mcell_index=GridT.bcell_start[grid_index] + id;
		const int imcell=GridT.which_bigcell[mcell_index];

		const int iat=GridT.which_atom[mcell_index];
		
		const int it=ucell.iat2it[iat];
		Atom* atom=&ucell.atoms[it];
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

		for(int ib=0; ib<pw.bxyz; ib++)
		{
			double *p=&psir_ylm[ib][block_index[id]];
			// meshcell_pos: z is the fastest
			const double dr[3] = {
				GridT.meshcell_pos[ib][0] + mt[0],
				GridT.meshcell_pos[ib][1] + mt[1],
				GridT.meshcell_pos[ib][2] + mt[2]};	

			double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);	// distance between atom and grid
			if(distance > (ORB.Phi[it].getRcut()- 1.0e-15)) 
			{
				cal_flag[ib][id]=false;
				ZEROS(p, block_size[id]);
			}
			else
			{
				cal_flag[ib][id]=true;
				
				std::vector<double> ylma;
				//if(distance[id] > GridT.orbital_rmax) continue;
				//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
				if (distance < 1.0E-9) distance += 1.0E-9;
				
				Ylm::sph_harm (	ucell.atoms[it].nwl,
						dr[0] / distance,
						dr[1] / distance,
						dr[2] / distance,
						ylma);
				// these parameters are about interpolation
				// because once we know the distance from atom to grid point,
				// we can get the parameters we need to do interpolation and
				// store them first!! these can save a lot of effort.
				const double position = distance / delta_r;
				const int ip = static_cast<int>(position);
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
						const Numerical_Orbital_Lm &philn = ORB.Phi[it].PhiLN(
								atom->iw2l[iw],
								atom->iw2n[iw]);
						phi=c1*philn.psi_uniform[ip]+c2*philn.dpsi_uniform[ip]
							+ c3*philn.psi_uniform[ip+1] + c4*philn.dpsi_uniform[ip+1];
					}
					*p=phi * ylma[atom->iw2_ylm[iw]];
				} // end iw
			}// end distance<=(ORB.Phi[it].getRcut()-1.0e-15)
		}// end ib
	}// end id
	return;
}

// can be done by GPU
void Gint_Gamma::cal_band_rho(
	const int na_grid, 
	const int LD_pool, 
	const int*const block_iw, 
	const int*const bsize, 
	const int*const colidx,
	const bool*const*const cal_flag, 
	const double*const*const psir_ylm, 
	double*const*const psir_DM,
	const int*const vindex)
{

    //parameters for dsymm, dgemm and ddot
    constexpr char side='L', uplo='U';
    constexpr char transa='N', transb='N';
    constexpr double alpha_symm=1, alpha_gemm=2, beta=1;    
    constexpr int inc=1;

    for(int is=0; is<NSPIN; ++is)
    {

        for (int ia1=0; ia1<na_grid; ++ia1)
        {
            const int iw1_lo=block_iw[ia1];

            //ia1==ia2, diagonal part
            // find the first ib and last ib for non-zeros cal_flag
            int first_ib=0, last_ib=0;
            for(int ib=0; ib<pw.bxyz; ++ib)
            {
                if(cal_flag[ib][ia1])
                {
                    first_ib=ib;
                    break;
                }
            }
            for(int ib=pw.bxyz-1; ib>=0; --ib)
            {
                if(cal_flag[ib][ia1])
                {
                    last_ib=ib+1;
                    break;
                }
            }
            const int ib_length=last_ib-first_ib;
            if(ib_length<=0) continue;

            int cal_num=0;
            for(int ib=first_ib; ib<last_ib; ++ib)
            {
                cal_num += cal_flag[ib][ia1];
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
                    if(cal_flag[ib][ia1])
                    {
                        dsymv_(&uplo, &bsize[ia1],
                            &alpha_symm, &LOC.DM[is][iw1_lo][iw1_lo], &GridT.lgd,
                            &psir_ylm[ib][colidx[ia1]], &inc,
                            &beta, &psir_DM[ib][colidx[ia1]], &inc);
                    }
                }
            }
            
            //OUT(ofs_running, "diagonal part of psir_DM done");
            for (int ia2=ia1+1; ia2<na_grid; ++ia2)
            {
                int first_ib=0, last_ib=0;
                for(int ib=0; ib<pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        first_ib=ib;
                        break;
                    }
                }
                for(int ib=pw.bxyz-1; ib>=0; --ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        last_ib=ib+1;
                        break;
                    }
                }
                const int ib_length=last_ib-first_ib;
                if(ib_length<=0) continue;

                int cal_pair_num=0;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
                }
                const int iw2_lo=block_iw[ia2];
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
                        if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
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
    
        double*const rhop = CHR.rho[is];
        for(int ib=0; ib<pw.bxyz; ++ib)
        {
            const double r = ddot_(&colidx[na_grid], psir_ylm[ib], &inc, psir_DM[ib], &inc);
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

	// I guess Peize add this, mohan 2021-01-31
	omp_init_lock(&lock);
    const double ne = this->gamma_charge();
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
			bool **cal_flag=new bool*[pw.bxyz];
			for(int i=0; i<pw.bxyz; ++i)
			{
				cal_flag[i]=new bool[max_size];
			}
		
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
		
						const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
		
						// get the value: how many atoms has orbital value on this grid.
						const int na_grid = GridT.how_many_atoms[ grid_index ];
						if(na_grid==0) continue;

						// it's a uniform grid to save orbital values, so the delta_r is a constant.
						const double delta_r = ORB.dr_uniform;						
						
						// here vindex refers to local potentials
						int* vindex = Gint_Gamma::get_vindex(ncyz, ibx, jby, kbz);	
						
						//------------------------------------------------------
						// band size: number of columns of a band
						//------------------------------------------------------------------
						int* block_size = get_bsize(na_grid, grid_index);
						
						//------------------------------------------------------
						// index of wave functions for each block
						//------------------------------------------------------
						int *block_iw = get_block_iw(na_grid, grid_index, this->max_size);

						int* block_index = get_colidx(na_grid, grid_index);

						// set up band matrix psir_ylm and psir_DM
						const int LD_pool=max_size*ucell.nwmax;
						Array_Pool psir_ylm(pw.bxyz, LD_pool);
						Array_Pool psir_DM(pw.bxyz, LD_pool);
						ZEROS(psir_DM.ptr_1D, pw.bxyz*LD_pool);
						
						this->cal_psir_ylm_rho(na_grid, grid_index, delta_r,
								block_index, block_size, 
								cal_flag, psir_ylm.ptr_2D);
						
						this->cal_band_rho(na_grid, LD_pool, block_iw, block_size, block_index,
								cal_flag, psir_ylm.ptr_2D, psir_DM.ptr_2D, vindex);

						free(vindex);			vindex=nullptr;
						free(block_size);		block_size=nullptr;
						free(block_iw);			block_iw=nullptr;
						free(block_index);		block_index=nullptr;
					}// k
				}// j
			}// i
			
			for(int i=0; i<pw.bxyz; i++)
			{
				delete[] cal_flag[i];
			}
			delete[] cal_flag;
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

    const double ne = sum * ucell.omega / pw.ncxyz;
    //xiaohui add 'OUT_LEVEL', 2015-09-16
    if(OUT_LEVEL != "m") OUT(ofs_running, "ne", ne);
    timer::tick("Gint_Gamma","gamma_charge",'I');

    return ne;
}

