#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"
#include "../src_lcao/global_fp.h" // mohan add 2021-01-30
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
	Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size, cal_flag);

    //evaluate psi and dpsi on grids
	Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::bigpw->bxyz, LD_pool);
	Gint_Tools::Array_Pool<double> dpsir_ylm_x(GlobalC::bigpw->bxyz, LD_pool);
	Gint_Tools::Array_Pool<double> dpsir_ylm_y(GlobalC::bigpw->bxyz, LD_pool);
	Gint_Tools::Array_Pool<double> dpsir_ylm_z(GlobalC::bigpw->bxyz, LD_pool);

	Gint_Tools::cal_dpsir_ylm(
		na_grid, grid_index, delta_r,
		block_index, block_size, 
		cal_flag,
		psir_ylm.ptr_2D,
		dpsir_ylm_x.ptr_2D,
		dpsir_ylm_y.ptr_2D,
		dpsir_ylm_z.ptr_2D
	);

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		Gint_Tools::Array_Pool<double> dpsix_DM(GlobalC::bigpw->bxyz, LD_pool);
		Gint_Tools::Array_Pool<double> dpsiy_DM(GlobalC::bigpw->bxyz, LD_pool);
		Gint_Tools::Array_Pool<double> dpsiz_DM(GlobalC::bigpw->bxyz, LD_pool);
		ModuleBase::GlobalFunc::ZEROS(dpsix_DM.ptr_1D, GlobalC::bigpw->bxyz*LD_pool);
		ModuleBase::GlobalFunc::ZEROS(dpsiy_DM.ptr_1D, GlobalC::bigpw->bxyz*LD_pool);
		ModuleBase::GlobalFunc::ZEROS(dpsiz_DM.ptr_1D, GlobalC::bigpw->bxyz*LD_pool);

		//calculating g_i,mu(r) = sum_nu rho_mu,nu d/dx_i psi_nu(r), x_i=x,y,z
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			Gint_Tools::mult_psi_DM(
				na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_x.ptr_2D,
				dpsix_DM.ptr_2D,
				inout->DM[is], 1);
			Gint_Tools::mult_psi_DM(
				na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_y.ptr_2D,
				dpsiy_DM.ptr_2D,
				inout->DM[is], 1);	
			Gint_Tools::mult_psi_DM(
				na_grid, LD_pool,
				block_iw, block_size,
				block_index, cal_flag,
				dpsir_ylm_z.ptr_2D,
				dpsiz_DM.ptr_2D,
				inout->DM[is], 1);
		}
		else
		{
			Gint_Tools::mult_psi_DMR(
				grid_index, na_grid,
				block_index, block_size,
				cal_flag, GlobalC::GridT,
				dpsir_ylm_x.ptr_2D,
				dpsix_DM.ptr_2D,
				inout->DM_R[is], 1);
			Gint_Tools::mult_psi_DMR(
				grid_index, na_grid,
				block_index, block_size,
				cal_flag, GlobalC::GridT,
				dpsir_ylm_y.ptr_2D,
				dpsiy_DM.ptr_2D,
				inout->DM_R[is], 1);
			Gint_Tools::mult_psi_DMR(
				grid_index, na_grid,
				block_index, block_size,
				cal_flag, GlobalC::GridT,
				dpsir_ylm_z.ptr_2D,
				dpsiz_DM.ptr_2D,
				inout->DM_R[is], 1);
		}

		//do sum_i,mu g_i,mu(r) * d/dx_i psi_mu(r) to get kinetic energy density on grid
		if(inout->job==Gint_Tools::job_type::tau)
		{
			this->cal_meshball_tau(
				na_grid, block_index,
				vindex,
				dpsir_ylm_x.ptr_2D, dpsir_ylm_y.ptr_2D, dpsir_ylm_z.ptr_2D,
				dpsix_DM.ptr_2D, dpsiy_DM.ptr_2D, dpsiz_DM.ptr_2D,
				inout->chr->kin_r[is]);
		}
	}

	delete[] block_iw;
	delete[] block_index;
	delete[] block_size;
	for(int ib=0; ib<GlobalC::bigpw->bxyz; ++ib)
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
	for(int ib=0; ib<GlobalC::bigpw->bxyz; ++ib)
	{
		double rx=ddot_(&block_index[na_grid], dpsix[ib], &inc, dpsix_dm[ib], &inc);
		double ry=ddot_(&block_index[na_grid], dpsiy[ib], &inc, dpsiy_dm[ib], &inc);
		double rz=ddot_(&block_index[na_grid], dpsiz[ib], &inc, dpsiz_dm[ib], &inc);
		const int grid = vindex[ib];
		rho[ grid ] += rx + ry + rz;
	}
}