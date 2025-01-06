#include "stress_func.h"
#include "module_base/math_integral.h"
#include "module_parameter/parameter.h"
#include "module_base/tool_threading.h"
#include "module_base/timer.h"
#include "module_base/libm/libm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

//calculate local pseudopotential stress in PW or VL_dVL stress in LCAO
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_loc(const UnitCell& ucell,
											 ModuleBase::matrix& sigma,
                                             ModulePW::PW_Basis* rho_basis,
											 const ModuleBase::matrix& vloc,
                                             const Structure_Factor* p_sf,
                                             const bool is_pw,
                                             const Charge* const chr)
{
    ModuleBase::TITLE("Stress_Func","stress_loc");
    ModuleBase::timer::tick("Stress_Func","stress_loc");

	std::vector<FPTYPE> dvloc(rho_basis->npw);
    FPTYPE evloc=0.0;
	FPTYPE fact=1.0;

	const int nspin_rho = (PARAM.inp.nspin == 2) ? 2 : 1;

	if (PARAM.inp.gamma_only && is_pw) 
	{ 
		fact=2.0;
	}

	std::vector<std::complex<FPTYPE>> aux(rho_basis->nmaxgr);

	/*
		blocking rho_basis->nrxx for data locality.

		By blocking aux with block size 1024,
		we can keep the blocked aux in L1 cache when iterating PARAM.inp.nspin loop
		performance will be better when number of atom is quite huge
	*/
	const int block_ir = 1024;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int irb = 0; irb < rho_basis->nrxx; irb += block_ir)
	{
		// calculate the actual task length of this block
		int ir_end = std::min(irb + block_ir, rho_basis->nrxx);

		{ // is = 0
			for (int ir = irb; ir < ir_end; ++ir)
			{ // initialize aux
				aux[ir] = std::complex<FPTYPE>(chr->rho[0][ir], 0.0 );
			}
		}
		if(nspin_rho == 2)
		{
			for (int ir = irb; ir < ir_end; ++ir)
			{ // accumulate aux
				aux[ir] += std::complex<FPTYPE>(chr->rho[1][ir], 0.0 );
			}
		}
 	}
	rho_basis->real2recip(aux.data(),aux.data());

//    if(INPUT.gamma_only==1) fact=2.0;
//    else fact=1.0;

	if(is_pw)
	{
#pragma omp parallel for collapse(2) reduction(+:evloc)
		for (int it=0; it<ucell.ntype; it++)
		{
			for (int ig=0; ig<rho_basis->npw; ig++)
			{
                if (rho_basis->ig_gge0 == ig) {
                    evloc += vloc(it, rho_basis->ig2igg[ig])
                             * (p_sf->strucFac(it, ig) * conj(aux[ig])).real();
                } else {
                    evloc += vloc(it, rho_basis->ig2igg[ig])
                             * (p_sf->strucFac(it, ig) * conj(aux[ig]) * fact).real();
}
            }
		}
    }
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const Atom* atom = &ucell.atoms[it];
		if(atom->coulomb_potential)
		{
		//
		// special case: pseudopotential is coulomb 1/r potential
		//
			this->dvloc_coulomb (ucell,atom->ncpp.zv, dvloc.data(), rho_basis);
		//
		}
		else
		{
		//
		// normal case: dvloc contains dV_loc(G)/dG
		//
			this->dvloc_of_g ( atom->ncpp.msh, atom->ncpp.rab.data(), atom->ncpp.r.data(),
					atom->ncpp.vloc_at.data(), atom->ncpp.zv, dvloc.data(), rho_basis, ucell);
		//
		}
#ifndef _OPENMP
		ModuleBase::matrix &local_sigma = sigma;
#else
#pragma omp parallel
{
		ModuleBase::matrix local_sigma(3, 3);
		#pragma omp for
#endif
		for(int ig = 0;ig< rho_basis->npw;ig++)
		{
			for (int l = 0;l< 3;l++)
			{
				for (int m = 0; m<l+1;m++)
				{
                    local_sigma(l, m) = local_sigma(l, m)
                                        + (conj(aux[ig]) * p_sf->strucFac(it, ig)).real() * 2.0
                                              * dvloc[rho_basis->ig2igg[ig]] * ucell.tpiba2
                                              * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] * fact;
                }
			}
		}
#ifdef _OPENMP
		#pragma omp critical(stress_loc_reduce)
		{
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<l+1;m++)
				{
					sigma(l,m) += local_sigma(l,m);
				}
			}
		}
}
#endif
	}

	if(is_pw)
	{
		for(int l = 0;l< 3;l++)
		{
			sigma(l,l) += evloc;
		}
	}
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
            Parallel_Reduce::reduce_pool(sigma(l, m));
		}
	}
	for(int l = 0;l< 3;l++)
	{
		for (int m = 0; m<l+1; m++)
		{
			sigma(m,l) = sigma(l,m);
		}
	}



	ModuleBase::timer::tick("Stress_Func","stress_loc");
	return;
}


template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::dvloc_of_g
(
const int& msh,
const FPTYPE* rab,
const FPTYPE* r,
const FPTYPE* vloc_at,
const FPTYPE& zp,
FPTYPE*  dvloc,
ModulePW::PW_Basis* rho_basis,
const UnitCell& ucell_in
)
{
  //----------------------------------------------------------------------
  //
  // dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  //

  //
  //FPTYPE  dvloc[ngl];
  // the fourier transform dVloc/dG
  //
	

	int igl0;
	this->device = base_device::get_device_type<Device>(this->ctx);

	std::vector<double> gx_arr(rho_basis->ngg+1);
    double* gx_arr_d = nullptr;
	// counter on erf functions or gaussians
	// counter on g shells vectors
	// first shell with g != 0
	
	std::vector<FPTYPE> aux(msh);

	// the  G=0 component is not computed
	if (rho_basis->gg_uniq[0] < 1.0e-8)
	{
		dvloc[0] = 0.0;
		igl0 = 1;
	}
	else
	{
		igl0 = 0;
	}

	// Pseudopotentials in numerical form (Vloc contains the local part)
	// In order to perform the Fourier transform, a term erf(r)/r is
	// subtracted in real space and added again in G space
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int igl = igl0; igl < rho_basis->ngg; igl++) {
        gx_arr[igl] = sqrt(rho_basis->gg_uniq[igl]) * ucell_in.tpiba;
    }

	gx_arr[rho_basis->ngg] = zp * ModuleBase::e2;

	//
	//   This is the part of the integrand function
	//   indipendent of |G| in real space
	//
	for(int i = 0;i< msh; i++)
	{
		aux[i] = r [i] * vloc_at [i] + zp * ModuleBase::e2 * erf(r[i]);
	}



    double *r_d = nullptr;
	double *rhoc_d = nullptr;
	double *rab_d = nullptr;
    double *aux_d = nullptr;
	double *drhocg_d = nullptr;
    if (this->device == base_device::GpuDevice) {
        resmem_var_op()(this->ctx, r_d, msh);
        resmem_var_op()(this->ctx, rhoc_d, msh);
        resmem_var_op()(this->ctx, rab_d, msh);

        resmem_var_op()(this->ctx, aux_d, msh);
        resmem_var_op()(this->ctx, gx_arr_d, rho_basis->ngg+1);
        resmem_var_op()(this->ctx, drhocg_d, rho_basis->ngg);

        syncmem_var_h2d_op()(this->ctx,
                             this->cpu_ctx,
                             gx_arr_d,
                             gx_arr.data(),
                             rho_basis->ngg+1);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, r_d, r, msh);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, rab_d, rab, msh);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, rhoc_d, aux.data(), msh);
    }



	if(this->device == base_device::GpuDevice) {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r_d,rhoc_d,gx_arr_d+igl0,rab_d,drhocg_d+igl0,msh,igl0,rho_basis->ngg-igl0,ucell_in.omega,3);
		syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, dvloc+igl0, drhocg_d+igl0, rho_basis->ngg-igl0);	

	} else {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r,aux.data(),gx_arr.data()+igl0,rab,dvloc+igl0,msh,igl0,rho_basis->ngg-igl0,ucell_in.omega,3);
	}	

	return;
}

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::dvloc_coulomb(const UnitCell& ucell, 
												const FPTYPE& zp, 
												FPTYPE* dvloc, 
												ModulePW::PW_Basis* rho_basis)
{
    int igl0;
	// start from |G|=0 or not.
    if (rho_basis->gg_uniq[0] < 1.0e-8)
    {
        dvloc[0] = 0.0;
        igl0 = 1;
    }
    else
    {
        igl0 = 0;
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = igl0; i < rho_basis->ngg; i++)
    {
        dvloc[i] = ModuleBase::FOUR_PI * zp * ModuleBase::e2 / ucell.omega
                   / pow((ucell.tpiba2 * rho_basis->gg_uniq[i]), 2);
    }

    return;
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif