#include "forces.h"
#include "stress_func.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"
// new
#include "module_base/complexmatrix.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/mathzone.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_general/module_vdw/vdw.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif



template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_cc(ModuleBase::matrix& forcecc,
                                          ModulePW::PW_Basis* rho_basis,
                                          const Charge* const chr)
{
    ModuleBase::TITLE("Forces", "cal_force_cc");
    // recalculate the exchange-correlation potential.
    ModuleBase::timer::tick("Forces", "cal_force_cc");

    int total_works = 0;
    // cal total works for skipping preprocess
    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            total_works += GlobalC::ucell.atoms[it].na;
        }
    }
    if (total_works == 0)
    {
        ModuleBase::timer::tick("Forces", "cal_force_cc");
        return;
    }

    ModuleBase::matrix v(GlobalV::NSPIN, rho_basis->nrxx);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
#ifdef USE_LIBXC
        const auto etxc_vtxc_v
            = XC_Functional::v_xc_meta(rho_basis->nrxx, GlobalC::ucell.omega, GlobalC::ucell.tpiba, chr);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc", "to use mGGA, compile with LIBXC");
#endif
    }
    else
    {
        if (GlobalV::NSPIN == 4)
            GlobalC::ucell.cal_ux();
        const auto etxc_vtxc_v = XC_Functional::v_xc(rho_basis->nrxx, chr, &GlobalC::ucell);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
    }

    const ModuleBase::matrix vxc = v;
    std::complex<double>* psiv = new std::complex<double>[rho_basis->nmaxgr];
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            psiv[ir] = std::complex<double>(vxc(0, ir), 0.0);
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            psiv[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
        }
    }

    // to G space
    rho_basis->real2recip(psiv, psiv);

    // psiv contains now Vxc(G)
    double* rhocg = new double[rho_basis->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocg, rho_basis->ngg);

    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {

            // chr->non_linear_core_correction(GlobalC::ppcell.numeric,
            //                                 GlobalC::ucell.atoms[it].ncpp.msh,
            //                                 GlobalC::ucell.atoms[it].ncpp.r,
            //                                 GlobalC::ucell.atoms[it].ncpp.rab,
            //                                 GlobalC::ucell.atoms[it].ncpp.rho_atc,
            //                                 rhocg);
            this->deriv_drhoc(GlobalC::ppcell.numeric,
                              GlobalC::ucell.atoms[it].ncpp.msh,
                              GlobalC::ucell.atoms[it].ncpp.r,
                              GlobalC::ucell.atoms[it].ncpp.rab,
                              GlobalC::ucell.atoms[it].ncpp.rho_atc,
                              rhocg,
                              rho_basis,
                              1);
#ifdef _OPENMP
#pragma omp parallel
            {
#endif
                for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ++ia)
                {
                    // get iat form table
                    int iat = GlobalC::ucell.itia2iat(it, ia);
                    double force[3] = {0, 0, 0};
#ifdef _OPENMP
#pragma omp for nowait
#endif
                    for (int ig = 0; ig < rho_basis->npw; ig++)
                    {
                        const ModuleBase::Vector3<double> gv = rho_basis->gcar[ig];
                        const ModuleBase::Vector3<double> pos = GlobalC::ucell.atoms[it].tau[ia];
                        const double rhocgigg = rhocg[rho_basis->ig2igg[ig]];
                        const std::complex<double> psiv_conj = conj(psiv[ig]);

                        const double arg = ModuleBase::TWO_PI * (gv.x * pos.x + gv.y * pos.y + gv.z * pos.z);
                        double sinp, cosp;
                        ModuleBase::libm::sincos(arg, &sinp, &cosp);
                        const std::complex<double> expiarg = std::complex<double>(sinp, cosp);

                        auto ipol0
                            = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.x * psiv_conj * expiarg;
                        force[0] += ipol0.real();

                        auto ipol1
                            = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.y * psiv_conj * expiarg;
                        force[1] += ipol1.real();

                        auto ipol2
                            = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.z * psiv_conj * expiarg;
                        force[2] += ipol2.real();
                    }
#ifdef _OPENMP
#pragma omp critical  
                    if (omp_get_num_threads() > 1)
                    {
                        forcecc(iat, 0) += force[0];
                        forcecc(iat, 1) += force[1];
                        forcecc(iat, 2) += force[2];
                    }
                    else
#endif
                    {
                        forcecc(iat, 0) += force[0];
                        forcecc(iat, 1) += force[1];
                        forcecc(iat, 2) += force[2];
                    }
                }
#ifdef _OPENMP
            } // omp parallel
#endif
        }
    }

    delete[] rhocg;

    delete[] psiv;                                                           // mohan fix bug 2012-03-22
    Parallel_Reduce::reduce_pool(forcecc.c, forcecc.nr * forcecc.nc); // qianrui fix a bug for kpar > 1
    ModuleBase::timer::tick("Forces", "cal_force_cc");
    return;
}



template<typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::deriv_drhoc
(
	const bool &numeric,
	const int mesh,
	const FPTYPE *r,
	const FPTYPE *rab,
	const FPTYPE *rhoc,
	FPTYPE *drhocg,
	ModulePW::PW_Basis* rho_basis,
	int type
)
{
	int  igl0;
	double gx = 0, rhocg1 = 0;
	//double *aux = new double[mesh];
	std::vector<double> aux(mesh);
	this->device = base_device::get_device_type<Device>(this->ctx);
	// the modulus of g for a given shell
	// the fourier transform
	// auxiliary memory for integration
	//double *gx_arr = new double[rho_basis->ngg];
	std::vector<double> gx_arr(rho_basis->ngg);
	double *gx_arr_d = nullptr;
	// counter on radial mesh points
	// counter on g shells
	// lower limit for loop on ngl

	//
	// G=0 term
	//
	if(type == 0){
		if (rho_basis->gg_uniq[0] < 1.0e-8)
		{
			drhocg [0] = 0.0;
			igl0 = 1;
		}
		else
		{
			igl0 = 0;
		}
	} else {
		if (rho_basis->gg_uniq[0] < 1.0e-8)
		{
			for (int ir = 0;ir < mesh; ir++)
			{
				aux [ir] = r [ir] * r [ir] * rhoc [ir];
			}
			ModuleBase::Integral::Simpson_Integral(mesh, aux.data(), rab, rhocg1);
			drhocg [0] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;
			igl0 = 1;
		} 
		else
		{
			igl0 = 0;
		}		
	}


	//
	// G <> 0 term
	//]

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int igl = igl0;igl< rho_basis->ngg;igl++)
	{
		gx_arr[igl] = sqrt(rho_basis->gg_uniq[igl] * GlobalC::ucell.tpiba2);
	}

	double *r_d = nullptr;
    double *rhoc_d = nullptr;
    double *rab_d = nullptr;
    double *aux_d = nullptr;
    double *drhocg_d = nullptr;
	if(this->device == base_device::GpuDevice ) {
		resmem_var_op()(this->ctx, r_d, mesh);
		resmem_var_op()(this->ctx, rhoc_d, mesh);
		resmem_var_op()(this->ctx, rab_d, mesh);

		resmem_var_op()(this->ctx, aux_d, mesh);
		resmem_var_op()(this->ctx, gx_arr_d, rho_basis->ngg);
		resmem_var_op()(this->ctx, drhocg_d, rho_basis->ngg);

		syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, gx_arr_d, gx_arr.data(), rho_basis->ngg);
		syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, r_d, r, mesh);
		syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, rab_d, rab, mesh);
		syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, rhoc_d, rhoc, mesh);
	}

	if(this->device == base_device::GpuDevice) {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r_d,rhoc_d,gx_arr_d+igl0,rab_d,drhocg_d+igl0,mesh,igl0,rho_basis->ngg-igl0,GlobalC::ucell.omega,type);
		syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, drhocg+igl0, drhocg_d+igl0, rho_basis->ngg-igl0);	



	} else {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r,rhoc,gx_arr.data()+igl0,rab,drhocg+igl0,mesh,igl0,rho_basis->ngg-igl0,GlobalC::ucell.omega,type);
    }

    delmem_var_op()(this->ctx, r_d);
    delmem_var_op()(this->ctx, rhoc_d);
    delmem_var_op()(this->ctx, rab_d);
    delmem_var_op()(this->ctx, gx_arr_d);
    delmem_var_op()(this->ctx, drhocg_d);
    return;
}


template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif