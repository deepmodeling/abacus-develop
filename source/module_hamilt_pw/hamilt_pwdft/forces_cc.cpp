#include "forces.h"
#include "stress_func.h"
#include "module_parameter/parameter.h"
#include "module_io/output_log.h"
// new
#include "module_base/complexmatrix.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/mathzone.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

#ifdef USE_LIBXC
#include "module_hamilt_general/module_xc/xc_functional_libxc.h"
#endif


template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_cc(ModuleBase::matrix& forcecc,
                                          ModulePW::PW_Basis* rho_basis,
                                          const Charge* const chr,
                                          const bool* numeric,
                                           UnitCell& ucell_in)
{
    ModuleBase::TITLE("Forces", "cal_force_cc");
    // recalculate the exchange-correlation potential.
    ModuleBase::timer::tick("Forces", "cal_force_cc");

    int total_works = 0;
    // cal total works for skipping preprocess
    for (int it = 0; it < ucell_in.ntype; ++it)
    {
        if (ucell_in.atoms[it].ncpp.nlcc)
        {
            total_works += ucell_in.atoms[it].na;
        }
    }
    if (total_works == 0)
    {
        ModuleBase::timer::tick("Forces", "cal_force_cc");
        return;
    }

    ModuleBase::matrix v(PARAM.inp.nspin, rho_basis->nrxx);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
#ifdef USE_LIBXC
        const auto etxc_vtxc_v
            = XC_Functional_Libxc::v_xc_meta(XC_Functional::get_func_id(), rho_basis->nrxx, ucell_in.omega, ucell_in.tpiba, chr);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc", "to use mGGA, compile with LIBXC");
#endif
    }
    else
    {
        elecstate::cal_ux(ucell_in);
        const auto etxc_vtxc_v = XC_Functional::v_xc(rho_basis->nrxx, chr, &ucell_in);

        // etxc = std::get<0>(etxc_vtxc_v);
        // vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
    }

    const ModuleBase::matrix vxc = v;
    std::complex<double>* psiv = new std::complex<double>[rho_basis->nmaxgr];
    if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
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

    std::vector<double> gv_x(rho_basis->npw);
    std::vector<double> gv_y(rho_basis->npw);
    std::vector<double> gv_z(rho_basis->npw);
    std::vector<double> rhocgigg_vec(rho_basis->npw);
    double *gv_x_d = nullptr;
    double *gv_y_d = nullptr;
    double *gv_z_d = nullptr;
    double *force_d = nullptr;
    double *rhocgigg_vec_d = nullptr;
    std::complex<FPTYPE>* psiv_d = nullptr;
    this->device = base_device::get_device_type<Device>(this->ctx);


#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        gv_x[ig] = rho_basis->gcar[ig].x;
        gv_y[ig] = rho_basis->gcar[ig].y;
        gv_z[ig] = rho_basis->gcar[ig].z;
    }

	if(this->device == base_device::GpuDevice ) {
		resmem_var_op()(gv_x_d, rho_basis->npw);
        resmem_var_op()(gv_y_d, rho_basis->npw);
        resmem_var_op()(gv_z_d, rho_basis->npw);
        resmem_var_op()(rhocgigg_vec_d, rho_basis->npw);
        resmem_complex_op()(psiv_d, rho_basis->nmaxgr);
        resmem_var_op()(force_d, 3);

		syncmem_var_h2d_op()(gv_x_d, gv_x.data(), rho_basis->npw);
        syncmem_var_h2d_op()(gv_y_d, gv_y.data(), rho_basis->npw);
        syncmem_var_h2d_op()(gv_z_d, gv_z.data(), rho_basis->npw);
        syncmem_complex_h2d_op()(psiv_d, psiv, rho_basis->nmaxgr);
	}


    for (int it = 0; it < ucell_in.ntype; ++it)
    {
        if (ucell_in.atoms[it].ncpp.nlcc)
        {

            // chr->non_linear_core_correction(numeric.numeric,
            //                                 ucell_in.atoms[it].ncpp.msh,
            //                                 ucell_in.atoms[it].ncpp.r,
            //                                 ucell_in.atoms[it].ncpp.rab,
            //                                 ucell_in.atoms[it].ncpp.rho_atc,
            //                                 rhocg);
            this->deriv_drhoc(numeric,
                              ucell_in.atoms[it].ncpp.msh,
                              ucell_in.atoms[it].ncpp.r.data(),
                              ucell_in.atoms[it].ncpp.rab.data(),
                              ucell_in.atoms[it].ncpp.rho_atc.data(),
                              rhocg,
                              rho_basis,
                              1,
                              ucell_in);

#ifdef _OPENMP
#pragma omp parallel for
#endif                              
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {
                rhocgigg_vec[ig] = rhocg[rho_basis->ig2igg[ig]];
            }

            if(this->device == base_device::GpuDevice ) {
                syncmem_var_h2d_op()(rhocgigg_vec_d, rhocgigg_vec.data(), rho_basis->npw);
            }
            for (int ia = 0; ia < ucell_in.atoms[it].na; ++ia)
            {
                const ModuleBase::Vector3<double> pos = ucell_in.atoms[it].tau[ia];
                // get iat form table
                int iat = ucell_in.itia2iat(it, ia);
                double force[3] = {0, 0, 0};

                if(this->device == base_device::GpuDevice ) {
                    syncmem_var_h2d_op()(force_d, force, 3);
                    hamilt::cal_force_npw_op<FPTYPE, Device>()(
                        psiv_d, gv_x_d, gv_y_d, gv_z_d, rhocgigg_vec_d, force_d, pos.x, pos.y, pos.z, 
                        rho_basis->npw, ucell_in.omega, ucell_in.tpiba
                    );      
                    syncmem_var_d2h_op()(force, force_d, 3);	          
                
                } else {
                    hamilt::cal_force_npw_op<FPTYPE, Device>()(
                        psiv, gv_x.data(), gv_y.data(), gv_z.data(), rhocgigg_vec.data(), force, pos.x, pos.y, pos.z, 
                        rho_basis->npw, ucell_in.omega, ucell_in.tpiba
                    );  
                }

                {
                    forcecc(iat, 0) += force[0];
                    forcecc(iat, 1) += force[1];
                    forcecc(iat, 2) += force[2];
                }
            }

        }
    }
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(gv_x_d);
        delmem_var_op()(gv_y_d);
        delmem_var_op()(gv_z_d);
        delmem_var_op()(force_d);
        delmem_var_op()(rhocgigg_vec_d);
        delmem_complex_op()(psiv_d);
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
	int type,
    const UnitCell& ucell_in
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
			drhocg [0] = ModuleBase::FOUR_PI * rhocg1 / ucell_in.omega;
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
		gx_arr[igl] = sqrt(rho_basis->gg_uniq[igl] * ucell_in.tpiba2);
	}

	double *r_d = nullptr;
    double *rhoc_d = nullptr;
    double *rab_d = nullptr;
    double *aux_d = nullptr;
    double *drhocg_d = nullptr;
	if(this->device == base_device::GpuDevice ) {
		resmem_var_op()(r_d, mesh);
		resmem_var_op()(rhoc_d, mesh);
		resmem_var_op()(rab_d, mesh);

		resmem_var_op()(aux_d, mesh);
		resmem_var_op()(gx_arr_d, rho_basis->ngg);
		resmem_var_op()(drhocg_d, rho_basis->ngg);

		syncmem_var_h2d_op()(gx_arr_d, gx_arr.data(), rho_basis->ngg);
		syncmem_var_h2d_op()(r_d, r, mesh);
		syncmem_var_h2d_op()(rab_d, rab, mesh);
		syncmem_var_h2d_op()(rhoc_d, rhoc, mesh);
	}

	if(this->device == base_device::GpuDevice) {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r_d,rhoc_d,gx_arr_d+igl0,rab_d,drhocg_d+igl0,mesh,igl0,rho_basis->ngg-igl0,ucell_in.omega,type);
		syncmem_var_d2h_op()(drhocg+igl0, drhocg_d+igl0, rho_basis->ngg-igl0);	



	} else {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r,rhoc,gx_arr.data()+igl0,rab,drhocg+igl0,mesh,igl0,rho_basis->ngg-igl0,ucell_in.omega,type);
    }

    delmem_var_op()(r_d);
    delmem_var_op()(rhoc_d);
    delmem_var_op()(rab_d);
    delmem_var_op()(gx_arr_d);
    delmem_var_op()(drhocg_d);
    return;
}


template class Forces<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, base_device::DEVICE_GPU>;
#endif