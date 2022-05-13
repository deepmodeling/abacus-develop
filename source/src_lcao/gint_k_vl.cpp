#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"
//#include <mkl_cblas.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

inline int find_offset(const int id1, const int id2, const int iat1, const int iat2,
				int* find_start, int* find_end)
{
	const int R1x=GlobalC::GridT.ucell_index2x[id1];
	const int R2x=GlobalC::GridT.ucell_index2x[id2];
	const int dRx=R1x-R2x;
	const int R1y=GlobalC::GridT.ucell_index2y[id1];
	const int R2y=GlobalC::GridT.ucell_index2y[id2];
	const int dRy=R1y-R2y;
	const int R1z=GlobalC::GridT.ucell_index2z[id1];
	const int R2z=GlobalC::GridT.ucell_index2z[id2];
	const int dRz=R1z-R2z;

	const int index=GlobalC::GridT.cal_RindexAtom(dRx, dRy, dRz, iat2);
	
	int offset=-1;
	for(int* find=find_start; find < find_end; ++find)
	{
		if( find[0] == index )
		{
			offset = find - find_start;
			break;
		}
	}

	assert(offset < GlobalC::GridT.nad[iat1]);
	return offset;
}

void Gint_k::cal_meshball_vlocal(
	int na_grid,
	int LD_pool,
	int grid_index, 
	int* block_size,
	int* block_index,
	int* block_iw,
	bool** cal_flag, 
	int* at, 
	double** psir_ylm,
	double** psir_vlbr3,
	double* pvpR)
{
	char transa='N', transb='T';
	double alpha=1, beta=1;
	int allnw=block_index[na_grid];

	int k=GlobalC::pw.bxyz;
	for(int ia1=0; ia1<na_grid; ++ia1)
	{
		//if(all_out_of_range[ia1]) continue;
		//const int iw1_lo=block_iw[ia1];
		const int idx1=block_index[ia1];
		int m=block_size[ia1];
		const int iat1=at[ia1];
		const int T1 = GlobalC::ucell.iat2it[iat1];
		const int mcell_index1 = GlobalC::GridT.bcell_start[grid_index] + ia1;
		const int id1 = GlobalC::GridT.which_unitcell[mcell_index1];
		const int DM_start = GlobalC::GridT.nlocstartg[iat1];
		// nad : how many adjacent atoms for atom 'iat'
		int* find_start = GlobalC::GridT.find_R2[iat1];
		int* find_end = GlobalC::GridT.find_R2[iat1] + GlobalC::GridT.nad[iat1];
		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int iat2=at[ia2];
			const int T2 = GlobalC::ucell.iat2it[iat2];
			if (iat1 <= iat2)
			{
    			int cal_num=0;
    			for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    			{
    				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
    				    ++cal_num;
    			}

    			if(cal_num==0) continue;
    			
                const int idx2=block_index[ia2];
        		int n=block_size[ia2];
				//const int I2 = GlobalC::ucell.iat2ia[iat2];
				const int mcell_index2 = GlobalC::GridT.bcell_start[grid_index] + ia2;
				const int id2 = GlobalC::GridT.which_unitcell[mcell_index2];
				int offset;
				offset=find_offset(id1, id2, iat1, iat2,
						find_start, find_end);

				const int iatw = DM_start + GlobalC::GridT.find_R2st[iat1][offset];	

			    if(cal_num>GlobalC::pw.bxyz/4)
			    {
					k=GlobalC::pw.bxyz;
					dgemm_(&transa, &transb, &n, &m, &k, &alpha,
						&psir_vlbr3[0][idx2], &LD_pool, 
						&psir_ylm[0][idx1], &LD_pool,
						&beta, &pvpR[iatw], &n);
				}
    			else
    			{
					for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
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

void Gint_k::cal_vlocal_k(const double *vrs1, const Grid_Technique &GridT, const int spin)
{
	ModuleBase::TITLE("Gint_k","cal_vlocal_k");

	if(!pvpR_alloc_flag)
	{
		ModuleBase::WARNING_QUIT("Gint_k::cal_vlocal_k","pvpR has not been allocated yet!");
	}
	else
	{
		ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[spin], GlobalC::GridT.nnrg);
	}

	ModuleBase::timer::tick("Gint_k","vlocal");
	const int max_size = GlobalC::GridT.max_atom;
	const double dv = GlobalC::ucell.omega/this->ncxyz;

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
#ifdef _OPENMP
			double* pvpR_reduced_thread;
        	pvpR_reduced_thread = new double[GlobalC::GridT.nnrg];
        	ModuleBase::GlobalFunc::ZEROS(pvpR_reduced_thread, GlobalC::GridT.nnrg);
#endif
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
				const int ibx=i*GlobalC::pw.bx;
				for(int j=0; j<nby; j++)
				{
					const int jby=j*GlobalC::pw.by;
					// count the z according to big box.
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
						
						//------------------------------------------------------------------
						// extract the local potentials.
						//------------------------------------------------------------------
						double *vldr3 = Gint_Tools::get_vldr3(vrs1, ncyz, ibx, jby, kbz, dv);

                        const Gint_Tools::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(
                                na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);

		#ifdef _OPENMP
						cal_meshball_vlocal(na_grid, LD_pool, grid_index, 
							block_size, block_index, block_iw, cal_flag, at,
							psir_ylm.ptr_2D, psir_vlbr3.ptr_2D, 
							pvpR_reduced_thread);
		#else
						cal_meshball_vlocal(na_grid, LD_pool, grid_index, 
							block_size, at, block_index, block_iw, cal_flag, at,
							psir_ylm.ptr_2D, psir_vlbr3.ptr_2D, 
							this->pvpR_reduced[spin]);
		#endif
						free(vldr3);		vldr3=nullptr;
                        delete[] block_iw;
                        delete[] block_index;
                        delete[] block_size;

						for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
							free(cal_flag[ib]);
						free(cal_flag);			cal_flag=nullptr;
					}// int k
				}// int j
			} // int i

#ifdef _OPENMP
			#pragma omp critical(cal_vl_k)
			for(int innrg=0; innrg<GlobalC::GridT.nnrg; innrg++)
			{
				pvpR_reduced[spin][innrg] += pvpR_reduced_thread[innrg];
			}
			delete[] pvpR_reduced_thread;
#endif
    	} // end omp
#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
	} // end of max_size

	ModuleBase::timer::tick("Gint_k","vlocal");
	return;
}