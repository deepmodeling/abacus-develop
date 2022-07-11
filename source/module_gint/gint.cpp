#include "gint.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"
#include "../src_pw/global.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

void Gint::cal_gint(Gint_inout *inout)
{
	ModuleBase::TITLE("Gint_interface","cal_gint");
    ModuleBase::timer::tick("Gint_interface", "cal_gint");

	const int max_size = GlobalC::GridT.max_atom;
	const int LD_pool = max_size*GlobalC::ucell.nwmax;
    const int lgd = GlobalC::GridT.lgd;

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
            //prepare some constants
			const int ncyz = GlobalC::rhopw->ny*GlobalC::rhopw->nplane; // mohan add 2012-03-25
			const double dv = GlobalC::ucell.omega/this->ncxyz;
			
			// it's a uniform grid to save orbital values, so the delta_r is a constant.
			const double delta_r = GlobalC::ORB.dr_uniform;

            if((inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta) && !GlobalV::GAMMA_ONLY_LOCAL)
            {
                if(!pvpR_alloc_flag)
                {
                    ModuleBase::WARNING_QUIT("Gint_interface::cal_gint","pvpR has not been allocated yet!");
                }
                else
                {
                    ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[inout->ispin], GlobalC::GridT.nnrg);
                }
            }

            //perpare auxiliary arrays to store thread-specific values
#ifdef _OPENMP
			double* pvpR_thread;
			if(inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta)
			{
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    pvpR_thread = new double[GlobalC::GridT.nnrg];
                    ModuleBase::GlobalFunc::ZEROS(pvpR_thread, GlobalC::GridT.nnrg);
                }
                if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
                {
                    pvpR_thread = new double[lgd*lgd];
                    ModuleBase::GlobalFunc::ZEROS(pvpR_thread, lgd*lgd);
                }
			}

			ModuleBase::matrix fvl_dphi_thread;
			ModuleBase::matrix svl_dphi_thread;
			if(inout->job==Gint_Tools::job_type::force)
			{
				if(inout->isforce)
				{
					fvl_dphi_thread.create(inout->fvl_dphi->nr,inout->fvl_dphi->nc);
					fvl_dphi_thread.zero_out();
				}
				if(inout->isstress)
				{
					svl_dphi_thread.create(inout->svl_dphi->nr,inout->svl_dphi->nc);
					svl_dphi_thread.zero_out();
				}
			}

    		#pragma omp for
#endif
            // entering the main loop of grid points
			for(int grid_index = 0; grid_index < GlobalC::bigpw->nbxx; grid_index++)
			{

				// get the value: how many atoms has orbital value on this grid.
				const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];

				if(na_grid==0) continue;				

				if(inout->job == Gint_Tools::job_type::rho)
				{
					//int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);
					int* vindex = Gint_Tools::get_vindex(GlobalC::GridT.start_ind[grid_index], ncyz);
					this->gint_kernel_rho(na_grid, grid_index, delta_r, vindex, LD_pool, inout);
					delete[] vindex;
				}
				else if(inout->job == Gint_Tools::job_type::tau || inout->job == Gint_Tools::job_type::crosstaus)
				{
					int* vindex = Gint_Tools::get_vindex(GlobalC::GridT.start_ind[grid_index], ncyz);
					this->gint_kernel_tau(na_grid, grid_index, delta_r, vindex, LD_pool, inout);
					delete[] vindex;					
				}
				else if(inout->job == Gint_Tools::job_type::force)
				{
					double* vldr3 = Gint_Tools::get_vldr3(inout->vl, GlobalC::GridT.start_ind[grid_index], ncyz, dv);
					double** DM_in;
					if(GlobalV::GAMMA_ONLY_LOCAL) DM_in = inout->DM[GlobalV::CURRENT_SPIN];
					if(!GlobalV::GAMMA_ONLY_LOCAL) DM_in = inout->DM_R;
					#ifdef _OPENMP
						this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool, 
							DM_in, inout->isforce, inout->isstress,
							&fvl_dphi_thread, &svl_dphi_thread);
					#else
						this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool, 
							DM_in, inout->isforce, inout->isstress,
							inout->fvl_dphi, inout->svl_dphi);
					#endif
					delete[] vldr3;
				}
				else if(inout->job==Gint_Tools::job_type::vlocal)
				{
					double* vldr3 = Gint_Tools::get_vldr3(inout->vl, GlobalC::GridT.start_ind[grid_index], ncyz, dv);
					#ifdef _OPENMP
						if((GlobalV::GAMMA_ONLY_LOCAL && lgd>0) || !GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool,
								pvpR_thread);
						}
					#else
						if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
						{
							this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool, pvpR_grid);
						}
						if(!GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool,
								this->pvpR_reduced[inout->ispin]);
						}
					#endif
					delete[] vldr3;
				}
				else if(inout->job==Gint_Tools::job_type::vlocal_meta)
				{
					double* vldr3 = Gint_Tools::get_vldr3(inout->vl, GlobalC::GridT.start_ind[grid_index], ncyz, dv);
					double* vkdr3 = Gint_Tools::get_vldr3(inout->vofk, GlobalC::GridT.start_ind[grid_index], ncyz, dv);
					#ifdef _OPENMP
						if((GlobalV::GAMMA_ONLY_LOCAL && lgd>0) || !GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool,
								pvpR_thread);
						}
					#else
						if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
						{
							this->gint_kernel_vlocal_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool, pvpR_grid);
						}
						if(!GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool,
								this->pvpR_reduced[inout->ispin]);
						}
					#endif
					delete[] vldr3;
					delete[] vkdr3;
				}				
			} // int grid_index

#ifdef _OPENMP
			if(inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta)
			{
                if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
                {
                    for(int i=0;i<lgd*lgd;i++)
                    {
                        #pragma omp critical(gint_gamma)
                        pvpR_grid[i] += pvpR_thread[i];
                    }
                    delete[] pvpR_thread;
                }
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    for(int innrg=0; innrg<GlobalC::GridT.nnrg; innrg++)
                    {
                        #pragma omp critical(gint_k)
                        pvpR_reduced[inout->ispin][innrg] += pvpR_thread[innrg];
                    }
                    delete[] pvpR_thread;
                }
			}

			if(inout->job==Gint_Tools::job_type::force)
			{
				if(inout->isforce)
				{
					#pragma omp critical(gint)
					inout->fvl_dphi[0]+=fvl_dphi_thread;
				}
				if(inout->isstress)
				{
					#pragma omp critical(gint)
					inout->svl_dphi[0]+=svl_dphi_thread;
				}
			}
#endif
		} // end of #pragma omp parallel

#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
    } // end of if (max_size)	

	ModuleBase::timer::tick("Gint_interface","cal_gint");
	return;
}

void Gint::prep_grid(
	const int &nbx_in,
	const int &nby_in,
	const int &nbz_in,
	const int &nbz_start_in,
    const int& ncxyz_in)
{
	ModuleBase::TITLE(GlobalV::ofs_running,"Gint_k","prep_grid");
    
    this->nbx = nbx_in;
	this->nby = nby_in;
	this->nbz = nbz_in;
	this->ncxyz = ncxyz_in;
	this->nbz_start = nbz_start_in;
	assert(nbx>0);
	assert(nby>0);
	assert(nbz>=0);
	assert(ncxyz>0);

	assert( GlobalC::ucell.omega > 0.0);

	return;
}