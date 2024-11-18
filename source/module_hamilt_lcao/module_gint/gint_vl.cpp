#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "gint_k.h"
#include "module_basis/module_ao/ORB_read.h"
#include "grid_technique.h"
#include "module_base/ylm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "module_base/array_pool.h"
//#include <mkl_cblas.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif


void Gint::cal_meshball_vlocal(
	const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int LD_pool,
	const int*const block_iw,				    // block_iw[na_grid],	index of wave functions for each block
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const int grid_index,                       // index of grid group, for tracing global atom index
	const bool*const*const cal_flag,	    	// cal_flag[this->bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	const double*const*const psir_ylm,		    // psir_ylm[this->bxyz][LD_pool]
	const double*const*const psir_vlbr3,	    // psir_vlbr3[this->bxyz][LD_pool]
	hamilt::HContainer<double>* hR)	    // this->hRGint is the container of <phi_0 | V | phi_R> matrix element.
{
	const char transa='N', transb='T';
	const double alpha=1, beta=1;
    const int lgd_now = this->gridt->lgd;

	const int mcell_index = this->gridt->bcell_start[grid_index];
	for(int ia1=0; ia1<na_grid; ++ia1)
	{
		const int bcell1 = mcell_index + ia1;
		const int iat1 = this->gridt->which_atom[bcell1];
		const int id1 = this->gridt->which_unitcell[bcell1];
		const int r1x = this->gridt->ucell_index2x[id1];
		const int r1y = this->gridt->ucell_index2y[id1];
		const int r1z = this->gridt->ucell_index2z[id1];

		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int bcell2 = mcell_index + ia2;
			const int iat2= this->gridt->which_atom[bcell2];
			const int id2 = this->gridt->which_unitcell[bcell2];
			const int r2x = this->gridt->ucell_index2x[id2];
			const int r2y = this->gridt->ucell_index2y[id2];
			const int r2z = this->gridt->ucell_index2z[id2];

			if(iat1<=iat2)
			{
                int first_ib=0;
                for(int ib=0; ib<this->bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        first_ib=ib;
                        break;
                    }
                }
                int last_ib=0;
                for(int ib=this->bxyz-1; ib>=0; --ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        last_ib=ib+1;
                        break;
                    }
                }
                const int ib_length = last_ib-first_ib;
                if(ib_length<=0) { continue; }

				// calculate the BaseMatrix of <iat1, iat2, R> atom-pair
				const int dRx = r1x - r2x;
            	const int dRy = r1y - r2y;
            	const int dRz = r1z - r2z;

				const auto tmp_matrix = hR->find_matrix(iat1, iat2, dRx, dRy, dRz);
				if (tmp_matrix == nullptr)
				{
					continue;
				}
				const int m = tmp_matrix->get_row_size();
				const int n = tmp_matrix->get_col_size();
                
				int cal_pair_num=0;
                for(int ib=first_ib;ib<last_ib; ++ib)
                {
                    cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
                }
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &n, &m, &ib_length, &alpha,
                        &psir_vlbr3[first_ib][block_index[ia2]], &LD_pool,
                        &psir_ylm[first_ib][block_index[ia1]], &LD_pool,
                        &beta, tmp_matrix->get_pointer(), &n); 
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
                                &beta, tmp_matrix->get_pointer(), &n);                          
                        }
                    }
                }
			}
		}
	}
}