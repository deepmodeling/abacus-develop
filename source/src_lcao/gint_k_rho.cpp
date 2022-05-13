#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"
#include "global_fp.h" // mohan add 2021-01-30

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

inline void cal_band_rho(
	const int size, 
	const int grid_index, 
	const int LD_pool, 
	int* block_iw, 
	int* block_size,
	int* block_index,
	int* at, 
	int* uc, 
	double** psir_ylm,
	int* vindex, 
    bool** cal_flag,
    double** DM_R,
	Charge* chr)
{
	char trans='N';
	double alpha_diag=1;
	double alpha_nondiag=2;
	double beta=1;
	int inc=1;
	int cal_num=0; // mohan initialize 2021-07-04
	int iw1_lo=0;

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		Gint_Tools::Array_Pool<double> psir_DM(GlobalC::pw.bxyz, LD_pool);
		ModuleBase::GlobalFunc::ZEROS(psir_DM.ptr_1D, GlobalC::pw.bxyz*LD_pool);

		for (int ia1=0; ia1<size; ++ia1)
		{
			const int iw1_lo=block_iw[ia1];
			const int iat1=at[ia1];
			const int id1=uc[ia1];
			const int idx1=block_index[ia1];
			const int R1x = GlobalC::GridT.ucell_index2x[id1];
			const int R1y = GlobalC::GridT.ucell_index2y[id1];
			const int R1z = GlobalC::GridT.ucell_index2z[id1];
			const int T1 = GlobalC::ucell.iat2it[iat1];
			int* find_start = GlobalC::GridT.find_R2[iat1];
			int* find_end = GlobalC::GridT.find_R2[iat1] + GlobalC::GridT.nad[iat1];
			//ia2==ia1
			cal_num=0;
			for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
			{
    			if(cal_flag[ib][ia1])
				{
    			    ++cal_num;
				}
			}

			int offset;
			if(cal_num>0)
			{
				//find offset				
				const int index = GlobalC::GridT.cal_RindexAtom(0, 0, 0, iat1);
				offset = -1;
				for(int* find=find_start; find < find_end; find++)
				{
					//--------------------------------------------------------------
					// start positions of adjacent atom of 'iat'
					//--------------------------------------------------------------
					if( find[0] == index ) 
					{
						offset = find - find_start; // start positions of adjacent atom of 'iat'
						break;
					}
				}

				assert(offset!=-1);
				assert(offset < GlobalC::GridT.nad[iat1]);				
			}

			if(cal_num>GlobalC::pw.bxyz/4)
			{				
				const int DM_start = GlobalC::GridT.nlocstartg[iat1]+ GlobalC::GridT.find_R2st[iat1][offset];					
				dgemm_(&trans, &trans, &block_size[ia1], &GlobalC::pw.bxyz, &block_size[ia1], &alpha_diag,
					&DM_R[is][DM_start], &block_size[ia1], 
					&psir_ylm[0][idx1], &LD_pool,  
					&beta, &psir_DM.ptr_2D[0][idx1], &LD_pool);
			}
			else if(cal_num>0)
			{	
				const int DM_start = GlobalC::GridT.nlocstartg[iat1]+ GlobalC::GridT.find_R2st[iat1][offset];
				for(int ib=0; ib<GlobalC::pw.bxyz; ++ib					)
				{
    				if(cal_flag[ib][ia1])
    				{
        				dgemv_(&trans, &block_size[ia1], &block_size[ia1], &alpha_diag,
					            &DM_R[is][DM_start], &block_size[ia1], 
					            &psir_ylm[ib][idx1], &inc,  
					            &beta, &psir_DM.ptr_2D[ib][idx1], &inc);
    				}
				}
			}

			//ia2>ia1
			for(int ia2=ia1+1; ia2<size; ++ia2)
			{			
			    cal_num=0;
    			for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    			{
        			if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        			    ++cal_num;
    			}
				
				int offset;
				if(cal_num>0)
				{
					const int iat2=at[ia2];
    				
    				// find offset
    				const int id2=uc[ia2];
    				const int R2x = GlobalC::GridT.ucell_index2x[id2];
    				const int R2y = GlobalC::GridT.ucell_index2y[id2];
    				const int R2z = GlobalC::GridT.ucell_index2z[id2];
    				const int dRx = R1x - R2x;
    				const int dRy = R1y - R2y;
    				const int dRz = R1z - R2z;
    				const int index = GlobalC::GridT.cal_RindexAtom(dRx, dRy, dRz, iat2);
    				offset = -1;
    				for(int* find=find_start; find < find_end; find++)
    				{
    					//--------------------------------------------------------------
    					// start positions of adjacent atom of 'iat'
    					//--------------------------------------------------------------
    					if( find[0] == index ) 
    					{
    						offset = find - find_start;
    						break;
    					}
    				}
					assert(offset!=-1);				
    				assert(offset < GlobalC::GridT.nad[iat1]);
				}

				if(cal_num>GlobalC::pw.bxyz/4)
				{
			        const int idx2=block_index[ia2];
    				const int DM_start = GlobalC::GridT.nlocstartg[iat1]+ GlobalC::GridT.find_R2st[iat1][offset];
    				dgemm_(&trans, &trans, &block_size[ia2], &GlobalC::pw.bxyz, &block_size[ia1], &alpha_nondiag,
    					&DM_R[is][DM_start], &block_size[ia2], 
    					&psir_ylm[0][idx1], &LD_pool,
    					&beta, &psir_DM.ptr_2D[0][idx2], &LD_pool);
				}
				else if(cal_num>0)
				{
					const int idx2=block_index[ia2];
    				const int DM_start = GlobalC::GridT.nlocstartg[iat1]+ GlobalC::GridT.find_R2st[iat1][offset];
					
    				for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    				{
        				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        				{
            				dgemv_(&trans, &block_size[ia2], &block_size[ia1], &alpha_nondiag,
            					&DM_R[is][DM_start], &block_size[ia2], 
            					&psir_ylm[ib][idx1], &inc,
            					&beta, &psir_DM.ptr_2D[ib][idx2], &inc);
        				}
    				}
				} // cal_num
			}// ia2
		} // ia1
		
		// calculate rho
		double *rhop = chr->rho[is];
		for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
		{
			double r=ddot_(&block_index[size], psir_ylm[ib], &inc, psir_DM.ptr_2D[ib], &inc);
			const int grid = vindex[ib];
			rhop[ grid ] += r;
		}
	}
}


void Gint_k::cal_rho_k(double** DM_R_in, Charge* chr)
{
	ModuleBase::TITLE("Gint_k","cal_rho_k");
    ModuleBase::timer::tick("Gint_k", "cal_rho_k");

    this->DM_R = DM_R_in;

	const int max_size = GlobalC::GridT.max_atom;

    if(max_size!=0)
    {
#ifdef __MKL
		const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(1);
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
			for(int i=0; i<nbx; i++)
			{
				const int ibx = i*GlobalC::pw.bx; // mohan add 2012-03-25
				for(int j=0; j<nby; j++)
				{
					const int jby = j*GlobalC::pw.by; // mohan add 2012-03-25
					for(int k=nbz_start; k<nbz_start+nbz; k++)
					{
						const int kbz = k*GlobalC::pw.bz-GlobalC::pw.nczp_start; //mohan add 2012-03-25
						
						const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;

						// get the value: how many atoms has orbital value on this grid.
						const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
						if(na_grid==0) continue;				
						
						// here vindex refers to local potentials
						int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);

                        int * block_iw, * block_index, * block_size, * at, * uc;
                        Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size, at, uc);

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

						cal_band_rho(
							na_grid, 
							grid_index, 
							LD_pool,
							block_iw, 
							block_size,
							block_index,
							at, 
							uc, 
							psir_ylm.ptr_2D,
							vindex, 
							cal_flag,
							DM_R,
							chr);
						
						free(vindex);			vindex=nullptr;
                        delete[] block_iw;
                        delete[] block_index;
                        delete[] block_size;

						for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
							free(cal_flag[ib]);
						free(cal_flag);			cal_flag=nullptr;
					}// int k
				}// int j
			} // int i
		} // end of #pragma omp parallel
#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
    } // end of if (max_size)	

	ModuleBase::timer::tick("Gint_k","cal_rho_k");
	return;
}


