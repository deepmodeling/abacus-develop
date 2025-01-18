#include "forces.h"
#include "module_io/output_log.h"
#include "stress_func.h"
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
void Forces<FPTYPE, Device>::cal_force_scc(ModuleBase::matrix& forcescc,
                                           ModulePW::PW_Basis* rho_basis,
                                           const ModuleBase::matrix& vnew,
                                           const bool vnew_exist,
                                           const bool* numeric,
                                           const UnitCell& ucell_in) {
    ModuleBase::TITLE("Forces", "cal_force_scc");
    ModuleBase::timer::tick("Forces", "cal_force_scc");

    // for orbital free case
    if (!vnew_exist) {
        ModuleBase::timer::tick("Forces", "cal_force_scc");
        return;
    }

    std::vector<std::complex<double>> psic(rho_basis->nmaxgr);

    const int nrxx = vnew.nc;
    const int nspin = vnew.nr;

    if (nspin == 1 || nspin == 4) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ir++) {
            psic[ir] = vnew(0, ir);
        }
    } else {
        int isup = 0;
        int isdw = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ir++) {
            psic[ir] = (vnew(isup, ir) + vnew(isdw, ir)) * 0.5;
        }
    }

    int ndm = 0;

    for (int it = 0; it < ucell_in.ntype; it++) {
        if (ndm < ucell_in.atoms[it].ncpp.msh) {
            ndm = ucell_in.atoms[it].ncpp.msh;
        }
    }

    // work space
	std::vector<double> rhocgnt(rho_basis->ngg);
    ModuleBase::GlobalFunc::ZEROS(rhocgnt.data(), rho_basis->ngg);

    rho_basis->real2recip(psic.data(), psic.data());

    int igg0 = 0;
    const int ig0 = rho_basis->ig_gge0;
    if (rho_basis->gg_uniq[0] < 1.0e-8) {
        igg0 = 1;
}

    double fact = 2.0;
    for (int nt = 0; nt < ucell_in.ntype; nt++) {
        //		Here we compute the G.ne.0 term
        const int mesh = ucell_in.atoms[nt].ncpp.msh;
        this->deriv_drhoc_scc(numeric,
                            mesh,
                            ucell_in.atoms[nt].ncpp.r.data(),
                            ucell_in.atoms[nt].ncpp.rab.data(),
                            ucell_in.atoms[nt].ncpp.rho_at.data(),
                            rhocgnt.data(),
                            rho_basis,
                            ucell_in);        
        int iat = 0;
        for (int it = 0; it < ucell_in.ntype; it++) {
            for (int ia = 0; ia < ucell_in.atoms[it].na; ia++) {
                if (nt == it) {
                    const ModuleBase::Vector3<double> pos
                        = ucell_in.atoms[it].tau[ia];
                    double &force0 = forcescc(iat, 0),
                           &force1 = forcescc(iat, 1),
                           &force2 = forcescc(iat, 2);
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : force0) reduction(+ : force1) reduction(+ : force2)
#endif
                    for (int ig = 0; ig < rho_basis->npw; ++ig) {
                        if (ig == ig0) {
                            continue;
}
                        const ModuleBase::Vector3<double> gv
                            = rho_basis->gcar[ig];
                        const double rhocgntigg
                            = rhocgnt[rho_basis->ig2igg[ig]];
                        const double arg = ModuleBase::TWO_PI * (gv * pos);
                        double sinp, cosp;
                        ModuleBase::libm::sincos(arg, &sinp, &cosp);
                        const std::complex<double> cpm
                            = std::complex<double>(sinp, cosp) * conj(psic[ig]);

                        force0 += fact * rhocgntigg * ucell_in.tpiba
                                  * gv.x * cpm.real();
                        force1 += fact * rhocgntigg * ucell_in.tpiba
                                  * gv.y * cpm.real();
                        force2 += fact * rhocgntigg * ucell_in.tpiba
                                  * gv.z * cpm.real();
                    }
                }
                iat++;
            }
        }
    }


    Parallel_Reduce::reduce_pool(forcescc.c, forcescc.nr * forcescc.nc);

    ModuleBase::timer::tick("Forces", "cal_force_scc");
    return;
}


template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::deriv_drhoc_scc(const bool& numeric,
                                              const int mesh,
                                              const FPTYPE* r,
                                              const FPTYPE* rab,
                                              const FPTYPE* rhoc,
                                              FPTYPE* drhocg,
                                              ModulePW::PW_Basis* rho_basis,
                                              const UnitCell& ucell_in) {
    int igl0 = 0;
    double gx = 0;
    double rhocg1 = 0;
    this->device = base_device::get_device_type<Device>(this->ctx);
    /// the modulus of g for a given shell
    /// the fourier transform
    /// auxiliary memory for integration
	std::vector<double> gx_arr(rho_basis->ngg);
    double* gx_arr_d = nullptr;
    /// counter on radial mesh points
    /// counter on g shells
    /// lower limit for loop on ngl

    ///
    /// G=0 term
    ///
    if (rho_basis->gg_uniq[0] < 1.0e-8) {
        drhocg[0] = 0.0;
        igl0 = 1;
    } else {
        igl0 = 0;
    }
    

    ///
    /// G <> 0 term
    ///]

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int igl = igl0; igl < rho_basis->ngg; igl++) {
        gx_arr[igl] = sqrt(rho_basis->gg_uniq[igl]) * ucell_in.tpiba;
    }

	double *r_d = nullptr;
    double *rhoc_d = nullptr;
    double *rab_d = nullptr;
    double *aux_d = nullptr;
    double *drhocg_d = nullptr;
    if (this->device == base_device::GpuDevice) {
        resmem_var_op()(r_d, mesh);
        resmem_var_op()(rhoc_d, mesh);
        resmem_var_op()(rab_d, mesh);

        resmem_var_op()(aux_d, mesh);
        resmem_var_op()(gx_arr_d, rho_basis->ngg);
        resmem_var_op()(drhocg_d, rho_basis->ngg);

        syncmem_var_h2d_op()(gx_arr_d,
                             gx_arr.data(),
                             rho_basis->ngg);
        syncmem_var_h2d_op()(r_d, r, mesh);
        syncmem_var_h2d_op()(rab_d, rab, mesh);
        syncmem_var_h2d_op()(rhoc_d, rhoc, mesh);
    }

	if(this->device == base_device::GpuDevice) {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r_d,rhoc_d,gx_arr_d+igl0,rab_d,drhocg_d+igl0,mesh,igl0,rho_basis->ngg-igl0,ucell_in.omega,2);
		syncmem_var_d2h_op()(drhocg+igl0, drhocg_d+igl0, rho_basis->ngg-igl0);	

	} else {
		hamilt::cal_stress_drhoc_aux_op<FPTYPE, Device>()(
			r,rhoc,gx_arr.data()+igl0,rab,drhocg+igl0,mesh,igl0,rho_basis->ngg-igl0,ucell_in.omega,2);

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