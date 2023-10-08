#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "gint_k.h"
#include "module_basis/module_ao/ORB_read.h"
#include "grid_technique.h"
#include "module_base/ylm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "gint_tools.h"

void Gint::gint_kernel_tau(
	const int na_grid,
	const int grid_index,
	const double delta_r,
	int* vindex,
	const int LD_pool,
	Gint_inout *inout)
{
	//prepare block information
	int * block_iw, * block_index, * block_size;
	bool** cal_flag;
	Gint_Tools::get_block_info(*this->gridt, this->bxyz, na_grid, grid_index, block_iw, block_index, block_size, cal_flag);

    //evaluate psi and dpsi on grids
	Gint_Tools::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
	Gint_Tools::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
	Gint_Tools::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
	Gint_Tools::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

	Gint_Tools::cal_dpsir_ylm(*this->gridt, 
		this->bxyz, na_grid, grid_index, delta_r,
		block_index, block_size, 
		cal_flag,
		psir_ylm.ptr_2D,
		dpsir_ylm_x.ptr_2D,
		dpsir_ylm_y.ptr_2D,
		dpsir_ylm_z.ptr_2D
	);

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		Gint_Tools::Array_Pool<double> dpsix_DM(this->bxyz, LD_pool);
		Gint_Tools::Array_Pool<double> dpsiy_DM(this->bxyz, LD_pool);
		Gint_Tools::Array_Pool<double> dpsiz_DM(this->bxyz, LD_pool);
		ModuleBase::GlobalFunc::ZEROS(dpsix_DM.ptr_1D, this->bxyz*LD_pool);
		ModuleBase::GlobalFunc::ZEROS(dpsiy_DM.ptr_1D, this->bxyz*LD_pool);
		ModuleBase::GlobalFunc::ZEROS(dpsiz_DM.ptr_1D, this->bxyz*LD_pool);

		//calculating g_i,mu(r) = sum_nu rho_mu,nu d/dx_i psi_nu(r), x_i=x,y,z
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			/*
			Gint_Tools::mult_psi_DM(
				*this->gridt,this->bxyz, na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_x.ptr_2D,
				dpsix_DM.ptr_2D,
				inout->DM[is], 1);
			Gint_Tools::mult_psi_DM(
				*this->gridt, this->bxyz, na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_y.ptr_2D,
				dpsiy_DM.ptr_2D,
				inout->DM[is], 1);	
			Gint_Tools::mult_psi_DM(
				*this->gridt, this->bxyz, na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_z.ptr_2D,
				dpsiz_DM.ptr_2D,
				inout->DM[is], 1);
			*/
			Gint_Tools::mult_psi_DM_new(
				*this->gridt,this->bxyz, grid_index, na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_x.ptr_2D,
				dpsix_DM.ptr_2D,
				this->DMRGint[is], 1);
			Gint_Tools::mult_psi_DM_new(
				*this->gridt, this->bxyz, grid_index, na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_y.ptr_2D,
				dpsiy_DM.ptr_2D,
				this->DMRGint[is], 1);	
			Gint_Tools::mult_psi_DM_new(
				*this->gridt, this->bxyz, grid_index, na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_z.ptr_2D,
				dpsiz_DM.ptr_2D,
				this->DMRGint[is], 1);
		}
		else
		{
			Gint_Tools::mult_psi_DMR(
				*this->gridt, this->bxyz, grid_index, na_grid,
				block_index, block_size,
				cal_flag, 
				dpsir_ylm_x.ptr_2D,
				dpsix_DM.ptr_2D,
				inout->DM_R[is],
				this->DMRGint[is],
				1);
			Gint_Tools::mult_psi_DMR(
				*this->gridt, this->bxyz, grid_index, na_grid,
				block_index, block_size,
				cal_flag,
				dpsir_ylm_y.ptr_2D,
				dpsiy_DM.ptr_2D,
				inout->DM_R[is], 
				this->DMRGint[is],
				1);
			Gint_Tools::mult_psi_DMR(
				*this->gridt, this->bxyz, grid_index, na_grid,
				block_index, block_size,
				cal_flag, 
				dpsir_ylm_z.ptr_2D,
				dpsiz_DM.ptr_2D,
				inout->DM_R[is], 
				this->DMRGint[is],
				1);
		}

		//do sum_i,mu g_i,mu(r) * d/dx_i psi_mu(r) to get kinetic energy density on grid
		if(inout->job==Gint_Tools::job_type::tau)
		{
			this->cal_meshball_tau(
				na_grid, block_index,
				vindex,
				dpsir_ylm_x.ptr_2D, dpsir_ylm_y.ptr_2D, dpsir_ylm_z.ptr_2D,
				dpsix_DM.ptr_2D, dpsiy_DM.ptr_2D, dpsiz_DM.ptr_2D,
				inout->rho[is]);
		}
	}

	delete[] block_iw;
	delete[] block_index;
	delete[] block_size;
	for(int ib=0; ib<this->bxyz; ++ib)
	{
		delete[] cal_flag[ib];
	}
	delete[] cal_flag;
}

void Gint::cal_meshball_tau(
	const int na_grid,
	int* block_index,
	int* vindex,
	double** dpsix,
	double** dpsiy,
	double** dpsiz,
	double** dpsix_dm,
	double** dpsiy_dm,
	double** dpsiz_dm,
	double* rho)
{		
	const int inc = 1;
	// sum over mu to get density on grid
	for(int ib=0; ib<this->bxyz; ++ib)
	{
		double rx=ddot_(&block_index[na_grid], dpsix[ib], &inc, dpsix_dm[ib], &inc);
		double ry=ddot_(&block_index[na_grid], dpsiy[ib], &inc, dpsiy_dm[ib], &inc);
		double rz=ddot_(&block_index[na_grid], dpsiz[ib], &inc, dpsiz_dm[ib], &inc);
		const int grid = vindex[ib];
		rho[ grid ] += rx + ry + rz;
	}
}