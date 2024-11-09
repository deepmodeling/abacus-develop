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


void Gint::cal_meshball_vlocal_gamma(
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
		const int iat1= this->gridt->which_atom[mcell_index + ia1];
		const int iw1_lo=block_iw[ia1];
		const int m=block_size[ia1];
		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int iat2= this->gridt->which_atom[mcell_index + ia2];
			const int iw2_lo=block_iw[ia2];
			if(iw1_lo<=iw2_lo)
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
				hamilt::AtomPair<double>* tmp_ap = hR->find_pair(iat1, iat2);
#ifdef __DEBUG
				assert(tmp_ap!=nullptr);
#endif
                int cal_pair_num=0;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
                }

                const int n=block_size[ia2];
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &n, &m, &ib_length, &alpha,
                        &psir_vlbr3[first_ib][block_index[ia2]], &LD_pool,
                        &psir_ylm[first_ib][block_index[ia1]], &LD_pool,
                        &beta, tmp_ap->get_pointer(0), &n);
						//&GridVlocal[iw1_lo*lgd_now+iw2_lo], &lgd_now);
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
                                &beta, tmp_ap->get_pointer(0), &n);
                        }
                    }
                }
			}
		}
	}
}

void Gint::cal_meshball_vlocal_k(
	const int na_grid,
	const int LD_pool,
	const int grid_index,
	const int*const block_size,
	const int*const block_index,
	const int*const block_iw,
	const bool*const*const cal_flag,
	const double*const*const psir_ylm,
	const double*const*const psir_vlbr3,
	double*const pvpR,
	const UnitCell &ucell)
{
    char transa = 'N', transb = 'T';
	double alpha=1, beta=1;
	int allnw=block_index[na_grid];

	int k=this->bxyz;
	for(int ia1=0; ia1<na_grid; ++ia1)
	{
		//if(all_out_of_range[ia1]) continue;
		//const int iw1_lo=block_iw[ia1];
		const int idx1=block_index[ia1];
		int m=block_size[ia1];
		const int mcell_index1 = this->gridt->bcell_start[grid_index] + ia1;
		const int iat1= this->gridt->which_atom[mcell_index1];
		const int T1 = ucell.iat2it[iat1];
		const int id1 = this->gridt->which_unitcell[mcell_index1];
		const int DM_start = this->gridt->nlocstartg[iat1];
		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int mcell_index2 = this->gridt->bcell_start[grid_index] + ia2;
			const int iat2 = this->gridt->which_atom[mcell_index2];
			const int T2 = ucell.iat2it[iat2];
			if (iat1 <= iat2)
			{
    			int cal_num=0;
    			for(int ib=0; ib<this->bxyz; ++ib)
    			{
    				if(cal_flag[ib][ia1] && cal_flag[ib][ia2]) {
    				    ++cal_num;
					}
    			}

    			if(cal_num==0) { continue; }

                const int idx2=block_index[ia2];
        		int n=block_size[ia2];
				//const int I2 = ucell.iat2ia[iat2];
				const int mcell_index2 = this->gridt->bcell_start[grid_index] + ia2;
				const int id2 = this->gridt->which_unitcell[mcell_index2];
				int offset;
				offset=this->gridt->find_offset(id1, id2, iat1, iat2);

				const int iatw = DM_start + this->gridt->find_R2st[iat1][offset];

			    if(cal_num>this->bxyz/4)
			    {
					k=this->bxyz;
					dgemm_(&transa, &transb, &n, &m, &k, &alpha,
						&psir_vlbr3[0][idx2], &LD_pool,
						&psir_ylm[0][idx1], &LD_pool,
						&beta, &pvpR[iatw], &n);
				}
    			else
    			{
					for(int ib=0; ib<this->bxyz; ++ib)
					{
						if(cal_flag[ib][ia1]&&cal_flag[ib][ia2])
						{
							k=1;
							dgemm_(&transa, &transb, &n, &m, &k, &alpha,
								&psir_vlbr3[ib][idx2], &LD_pool,
								&psir_ylm[ib][idx1], &LD_pool,
								&beta, &pvpR[iatw], &n);
						}
					}
    			}
			}
		}
	}
}