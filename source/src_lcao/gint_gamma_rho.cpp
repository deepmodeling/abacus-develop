//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

// can be done by GPU
void Gint_Gamma::cal_band_rho(
	const int na_grid,    							// how many atoms on this (i,j,k) grid
	const int LD_pool, 
    const int*const block_iw, 						// block_iw[na_grid],	index of wave functions for each block
    const int*const block_size, 					// block_size[na_grid],	band size: number of columns of a band
    const int*const block_index,					// block_index[na_grid+1], count total number of atomis orbitals
    const bool*const*const cal_flag, 				// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    const double*const*const psir_ylm,				// psir_ylm[GlobalC::pw.bxyz][LD_pool]
    const int*const vindex,							// vindex[GlobalC::pw.bxyz]
    const double*const*const*const DM,				// DM[GlobalV::NSPIN][lgd_now][lgd_now]
    Charge* chr) const		// rho[GlobalV::NSPIN][GlobalC::pw.nrxx]
{
    //parameters for dsymm, dgemm and ddot
    constexpr char side='L', uplo='U';
    constexpr char transa='N', transb='N';
    constexpr double alpha_symm=1, alpha_gemm=2, beta=1;    
    constexpr int inc=1;

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        Gint_Tools::Array_Pool<double> psir_DM(GlobalC::pw.bxyz, LD_pool);
        ModuleBase::GlobalFunc::ZEROS(psir_DM.ptr_1D, GlobalC::pw.bxyz*LD_pool);

        for (int ia1=0; ia1<na_grid; ++ia1)
        {
            const int iw1_lo=block_iw[ia1];

            //ia1==ia2, diagonal part
            // find the first ib and last ib for non-zeros cal_flag
            int first_ib=0, last_ib=0;
            for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
            {
                if(cal_flag[ib][ia1])
                {
                    first_ib=ib;
                    break;
                }
            }
            for(int ib=GlobalC::pw.bxyz-1; ib>=0; --ib)
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
                    &alpha_symm, &DM[is][iw1_lo][iw1_lo], &GlobalC::GridT.lgd, 
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
                            &alpha_symm, &DM[is][iw1_lo][iw1_lo], &GlobalC::GridT.lgd,
                            &psir_ylm[ib][block_index[ia1]], &inc,
                            &beta, &psir_DM.ptr_2D[ib][block_index[ia1]], &inc);
                    }
                }
            }
            
            //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "diagonal part of psir_DM done");
            for (int ia2=ia1+1; ia2<na_grid; ++ia2)
            {
                int first_ib=0, last_ib=0;
                for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        first_ib=ib;
                        break;
                    }
                }
                for(int ib=GlobalC::pw.bxyz-1; ib>=0; --ib)
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
                        &alpha_gemm, &DM[is][iw1_lo][iw2_lo], &GlobalC::GridT.lgd, 
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
                                &alpha_gemm, &DM[is][iw1_lo][iw2_lo], &GlobalC::GridT.lgd,
                                &psir_ylm[ib][block_index[ia1]], &inc,
                                &beta, &psir_DM.ptr_2D[ib][block_index[ia2]], &inc);
                        }
                    }
                }
                //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "upper triangle part of psir_DM done, atom2", ia2);
            }// ia2
        } // ia1
    
        double *rhop = chr->rho[is];
        for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
        {
            const double r = ddot_(&block_index[na_grid], psir_ylm[ib], &inc, psir_DM.ptr_2D[ib], &inc);
            const int grid = vindex[ib];
            rhop[grid] += r;
        }
    } // end is
}

// calculate charge density
void Gint_Gamma::cal_rho(double*** DM_in, Charge* chr)
{
    ModuleBase::TITLE("Gint_Gamma","cal_rho");
    ModuleBase::timer::tick("Gint_Gamma","cal_rho");

    max_size = GlobalC::GridT.max_atom;

	if(max_size)
    {
#ifdef __MKL
   		const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(std::max(1,mkl_threads/GlobalC::GridT.nbx));		// Peize Lin update 2021.01.20
#endif
		
#ifdef _OPENMP
		#pragma omp parallel
#endif
		{		
			const int nbx = GlobalC::GridT.nbx;
			const int nby = GlobalC::GridT.nby;
			const int nbz_start = GlobalC::GridT.nbzp_start;
			const int nbz = GlobalC::GridT.nbzp;
		
			const int ncyz = GlobalC::pw.ncy*GlobalC::pw.nczp; // mohan add 2012-03-25
            
            // it's a uniform grid to save orbital values, so the delta_r is a constant.
            const double delta_r = GlobalC::ORB.dr_uniform;		
#ifdef _OPENMP
			#pragma omp for
#endif
			for (int i=0; i<nbx; i++)
			{
				const int ibx = i*GlobalC::pw.bx;
				for (int j=0; j<nby; j++)
				{
					const int jby = j*GlobalC::pw.by;
					for (int k=nbz_start; k<nbz_start+nbz; k++)
					{
						const int kbz = k*GlobalC::pw.bz-GlobalC::pw.nczp_start;
		
						const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
		
						// get the value: how many atoms has orbital value on this grid.
						const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
						if(na_grid==0) continue;				
						
						// here vindex refers to local potentials
						int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);	
						
                        int * block_iw, * block_index, * block_size;
                        Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);

						//------------------------------------------------------
						// whether the atom-grid distance is larger than cutoff
						//------------------------------------------------------
						bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

						// set up band matrix psir_ylm and psir_DM
						const int LD_pool = max_size*GlobalC::ucell.nwmax;
						
						Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::cal_psir_ylm(
							na_grid, grid_index, delta_r,
							block_index, block_size, 
							cal_flag,
                            psir_ylm.ptr_2D);
						
						this->cal_band_rho(na_grid, LD_pool, block_iw, block_size, block_index,
							cal_flag, psir_ylm.ptr_2D, vindex, DM_in, chr);

						free(vindex);			vindex=nullptr;
                        delete[] block_iw;
                        delete[] block_index;
                        delete[] block_size;

						for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
							free(cal_flag[ib]);
						free(cal_flag);			cal_flag=nullptr;
					}// k
				}// j
			}// i
		} // end of #pragma omp parallel
			
#ifdef __MKL
   		mkl_set_num_threads(mkl_threads);
#endif
    } // end of if(max_size)

    ModuleBase::timer::tick("Gint_Gamma","cal_rho");
    return;
}
