#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "module_ORB/ORB_read.h"
#include "src_pw/global.h"
#include "module_base/blas_connector.h"

#include "global_fp.h" // mohan add 2021-01-30

// can be done by GPU
void Gint_Gamma::cal_band_rho(
	const int na_grid, 
	const int LD_pool, 
	const int*const block_iw, 
	const int*const block_size, 
	const int*const block_index,
	const bool*const*const cal_flag, 
	const double*const*const psir_ylm,
	const int*const vindex)
{
    //parameters for dsymm, dgemm and ddot
    constexpr char side='L', uplo='U';
    constexpr char transa='N', transb='N';
    constexpr double alpha_symm=1, alpha_gemm=2, beta=1;    
    constexpr int inc=1;

    for(int is=0; is<NSPIN; ++is)
    {
        Gint_Tools::Array_Pool<double> psir_DM(pw.bxyz, LD_pool);
        ZEROS(psir_DM.ptr_1D, pw.bxyz*LD_pool);

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
                dsymm_(&side, &uplo, &block_size[ia1], &ib_length, 
                    &alpha_symm, &LOC.DM[is][iw1_lo][iw1_lo], &GridT.lgd, 
                    &psir_ylm[first_ib][block_index[ia1]], &LD_pool, 
                    &beta, &psir_DM.ptr_2D[first_ib][block_index[ia1]], &LD_pool);
            }
            else
            {
                // int k=1;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    if(cal_flag[ib][ia1])
                    {
                        dsymv_(&uplo, &block_size[ia1],
                            &alpha_symm, &LOC.DM[is][iw1_lo][iw1_lo], &GridT.lgd,
                            &psir_ylm[ib][block_index[ia1]], &inc,
                            &beta, &psir_DM.ptr_2D[ib][block_index[ia1]], &inc);
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
                    dgemm_(&transa, &transb, &block_size[ia2], &ib_length, &block_size[ia1], 
                        &alpha_gemm, &LOC.DM[is][iw1_lo][iw2_lo], &GridT.lgd, 
                        &psir_ylm[first_ib][block_index[ia1]], &LD_pool, 
                        &beta, &psir_DM.ptr_2D[first_ib][block_index[ia2]], &LD_pool);
                }
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        {
                            dgemv_(&transa, &block_size[ia2], &block_size[ia1], 
                                &alpha_gemm, &LOC.DM[is][iw1_lo][iw2_lo], &GridT.lgd,
                                &psir_ylm[ib][block_index[ia1]], &inc,
                                &beta, &psir_DM.ptr_2D[ib][block_index[ia2]], &inc);
                        }
                    }
                }
                //OUT(ofs_running, "upper triangle part of psir_DM done, atom2", ia2);
            }// ia2
        } // ia1
    
        double*const rhop = CHR.rho[is];
        for(int ib=0; ib<pw.bxyz; ++ib)
        {
            const double r = ddot_(&block_index[na_grid], psir_ylm[ib], &inc, psir_DM.ptr_2D[ib], &inc);
            const int grid = vindex[ib];
            rhop[ grid ] += r;
        }
    } // end is
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
        const int omp_threads = omp_get_max_threads();
		omp_set_num_threads(std::max(1,omp_threads/GridT.nbx));			// Peize Lin update 2021.01.20
		
#ifdef __OPENMP
		#pragma omp parallel
#endif
		{		
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
						int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);	
						
						//------------------------------------------------------
						// band size: number of columns of a band
						//------------------------------------------------------------------
						int* block_size = Gint_Tools::get_block_size(na_grid, grid_index);
						
						//------------------------------------------------------
						// index of wave functions for each block
						//------------------------------------------------------
						int *block_iw = Gint_Tools::get_block_iw(na_grid, grid_index, this->max_size);

						int* block_index = Gint_Tools::get_block_index(na_grid, grid_index);

						//------------------------------------------------------
						// whether the atom-grid distance is larger than cutoff
						//------------------------------------------------------
						bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

						// set up band matrix psir_ylm and psir_DM
						const int LD_pool = max_size*ucell.nwmax;
						
						const Gint_Tools::Array_Pool<double> psir_ylm = Gint_Tools::cal_psir_ylm(
							na_grid, LD_pool, grid_index, delta_r,
							block_index, block_size, 
							cal_flag);
						
						this->cal_band_rho(na_grid, LD_pool, block_iw, block_size, block_index,
							cal_flag, psir_ylm.ptr_2D, vindex);

						free(vindex);			vindex=nullptr;
						free(block_size);		block_size=nullptr;
						free(block_iw);			block_iw=nullptr;
						free(block_index);		block_index=nullptr;

						for(int ib=0; ib<pw.bxyz; ++ib)
							free(cal_flag[ib]);
						free(cal_flag);			cal_flag=nullptr;
					}// k
				}// j
			}// i
		} // end of #pragma omp parallel
		
        for(int is=0; is<NSPIN; is++)
		{
            for (int ir=0; ir<pw.nrxx; ir++)
			{
                sum += CHR.rho[is][ir];
			}
		}
			
        omp_set_num_threads(omp_threads);
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

