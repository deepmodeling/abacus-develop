#include "stress_func.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/libm/libm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "kernels/stress_op.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//calcualte the Ewald stress term in PW and LCAO
template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_ewa(const UnitCell& ucell,
											 ModuleBase::matrix& sigma, 
											 ModulePW::PW_Basis* rho_basis, 
											 const bool is_pw)
{
    ModuleBase::TITLE("Stress_Func","stress_ewa");
    ModuleBase::timer::tick("Stress_Func","stress_ewa");

    FPTYPE charge=0;
    for(int it=0; it < ucell.ntype; it++)
	{
		charge = charge + ucell.atoms[it].ncpp.zv * ucell.atoms[it].na;
	}
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    FPTYPE alpha=2.9;
    FPTYPE upperbound;
    do{
       alpha-=0.1;
       if(alpha==0.0)
          ModuleBase::WARNING_QUIT("stres_ew", "optimal alpha not found");
       upperbound =ModuleBase::e2 * pow(charge,2) * sqrt( 2 * alpha / (ModuleBase::TWO_PI)) * erfc(sqrt(ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    }
    while(upperbound>1e-7);

    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    FPTYPE sdewald;
	const int ig0 = rho_basis->ig_gge0;
    if( ig0 >= 0)
	{
       sdewald = (ModuleBase::TWO_PI) * ModuleBase::e2 / 4.0 / alpha * pow(charge/ucell.omega,2);
    }
    else 
	{
       sdewald = 0.0;
    }

    //sdewald is the diagonal term 

    FPTYPE fact=1.0;
    if (PARAM.globalv.gamma_only_pw && is_pw) fact=2.0;
//    else fact=1.0;

    // Prepare data for the sincos op
    std::vector<FPTYPE> zv_facts_host(ucell.nat);
    std::vector<FPTYPE> tau_flat(ucell.nat * 3);

    for (int iat = 0; iat < ucell.nat; iat++) {
        int it = ucell.iat2it[iat];
        int ia = ucell.iat2ia[iat];

        zv_facts_host[iat] = static_cast<FPTYPE>(ucell.atoms[it].ncpp.zv);

        tau_flat[iat * 3 + 0] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia][0]);
        tau_flat[iat * 3 + 1] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia][1]);
        tau_flat[iat * 3 + 2] = static_cast<FPTYPE>(ucell.atoms[it].tau[ia][2]);
    }

    std::vector<FPTYPE> gcar_flat(rho_basis->npw * 3);
    for (int ig = 0; ig < rho_basis->npw; ig++) {
        gcar_flat[ig * 3 + 0] = static_cast<FPTYPE>(rho_basis->gcar[ig][0]);
        gcar_flat[ig * 3 + 1] = static_cast<FPTYPE>(rho_basis->gcar[ig][1]);
        gcar_flat[ig * 3 + 2] = static_cast<FPTYPE>(rho_basis->gcar[ig][2]);
    }

    // Allocate result arrays
    std::vector<FPTYPE> rhostar_real_host(rho_basis->npw);
    std::vector<FPTYPE> rhostar_imag_host(rho_basis->npw);

    // Device pointers for GPU data transfer
    FPTYPE* d_gcar = nullptr;
    FPTYPE* d_tau = nullptr;
    FPTYPE* d_zv_facts = nullptr;
    FPTYPE* d_rhostar_real = nullptr;
    FPTYPE* d_rhostar_imag = nullptr;

    // GPU memory management and data transfer
    if (this->device == base_device::GpuDevice) {
        // Allocate GPU memory
        resmem_var_op()(this->ctx, d_gcar, rho_basis->npw * 3);
        resmem_var_op()(this->ctx, d_tau, ucell.nat * 3);
        resmem_var_op()(this->ctx, d_zv_facts, ucell.nat);
        resmem_var_op()(this->ctx, d_rhostar_real, rho_basis->npw);
        resmem_var_op()(this->ctx, d_rhostar_imag, rho_basis->npw);

        // Data transfer H2D
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_gcar, gcar_flat.data(), rho_basis->npw * 3);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_tau, tau_flat.data(), ucell.nat * 3);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_zv_facts, zv_facts_host.data(), ucell.nat);

        // Initialize output arrays to zero on GPU
        base_device::memory::set_memory_op<FPTYPE, Device>()(this->ctx, d_rhostar_real, 0.0, rho_basis->npw);
        base_device::memory::set_memory_op<FPTYPE, Device>()(this->ctx, d_rhostar_imag, 0.0, rho_basis->npw);
    } else {
        // CPU case: use host pointers directly and initialize arrays to zero
        d_gcar = gcar_flat.data();
        d_tau = tau_flat.data();
        d_zv_facts = zv_facts_host.data();
        d_rhostar_real = rhostar_real_host.data();
        d_rhostar_imag = rhostar_imag_host.data();

        // Initialize output arrays to zero for CPU case
        std::fill(rhostar_real_host.begin(), rhostar_real_host.end(), static_cast<FPTYPE>(0.0));
        std::fill(rhostar_imag_host.begin(), rhostar_imag_host.end(), static_cast<FPTYPE>(0.0));
    }

    // Call sincos op (outside OpenMP parallel region, op has its own parallelization)
    hamilt::cal_stress_ewa_sincos_op<FPTYPE, Device>()(
        this->ctx,
        ucell.nat,
        rho_basis->npw,
        rho_basis->ig_gge0,
        d_gcar,
        d_tau,
        d_zv_facts,
        d_rhostar_real,
        d_rhostar_imag
    );

    // GPU data transfer D2H and memory cleanup
    if (this->device == base_device::GpuDevice) {
        // Transfer results back to host
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, rhostar_real_host.data(), d_rhostar_real, rho_basis->npw);
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, rhostar_imag_host.data(), d_rhostar_imag, rho_basis->npw);

        // Free GPU memory
        delmem_var_op()(this->ctx, d_gcar);
        delmem_var_op()(this->ctx, d_tau);
        delmem_var_op()(this->ctx, d_zv_facts);
        delmem_var_op()(this->ctx, d_rhostar_real);
        delmem_var_op()(this->ctx, d_rhostar_imag);
    }

#ifdef _OPENMP
#pragma omp parallel
{
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();
	ModuleBase::matrix local_sigma(3, 3);
	FPTYPE local_sdewald = 0.0;
#else
    int num_threads = 1;
    int thread_id = 0;
	ModuleBase::matrix& local_sigma = sigma;
	FPTYPE& local_sdewald = sdewald;
#endif

	// Calculate ig range of this thread, avoid thread sync
	int ig, ig_end;
	ModuleBase::TASK_DIST_1D(num_threads, thread_id, rho_basis->npw, ig, ig_end);
	ig_end = ig + ig_end;

    FPTYPE g2,g2a;
    FPTYPE sewald;
    for(; ig < ig_end; ig++)
	{
		if(ig == ig0)  continue;
		g2 = rho_basis->gg[ig]* ucell.tpiba2;
		g2a = g2 /4.0/alpha;

		// Use precomputed rhostar values
		FPTYPE rhostar_real = rhostar_real_host[ig] / ucell.omega;
		FPTYPE rhostar_imag = rhostar_imag_host[ig] / ucell.omega;

		// Calculate |rhostar|² - mathematically equivalent to pow(std::abs(rhostar), 2)
		FPTYPE rhostar_abs2 = rhostar_real * rhostar_real + rhostar_imag * rhostar_imag;

		sewald = fact* (ModuleBase::TWO_PI) * ModuleBase::e2 * ModuleBase::libm::exp(-g2a) / g2 * rhostar_abs2;
		local_sdewald -= sewald;
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				local_sigma(l, m) += sewald * ucell.tpiba2 * 2.0 * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / g2 * (g2a + 1);
			}
		}
	}

    //R-space sum here (only for the processor that contains G=0) 
    int mxr = 200;
    int *irr;
    ModuleBase::Vector3<FPTYPE> *r;
    FPTYPE *r2;
    FPTYPE rr;
    ModuleBase::Vector3<FPTYPE> d_tau;
    FPTYPE r0[3];
    FPTYPE rmax=0.0;
    int nrm=0;
    FPTYPE fac;
	if(ig0 >= 0)
	{
		r = new ModuleBase::Vector3<FPTYPE>[mxr];
		r2 = new FPTYPE[mxr];
		irr = new int[mxr];

		FPTYPE sqa = sqrt(alpha);
		FPTYPE sq8a_2pi = sqrt(8 * alpha / (ModuleBase::TWO_PI));
		rmax = 4.0/sqa/ucell.lat0;

		// collapse it, ia, jt, ja loop into a single loop
		long long ijat, ijat_end;
		int it, i, jt, j;
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, (long long)ucell.nat * ucell.nat, ijat, ijat_end);
		ijat_end = ijat + ijat_end;
		ucell.ijat2iaitjajt(ijat, &i, &it, &j, &jt);

		while (ijat < ijat_end)
		{
			if (ucell.atoms[it].na != 0 && ucell.atoms[jt].na != 0)
			{
				//calculate tau[na]-tau[nb]
				d_tau = ucell.atoms[it].tau[i] - ucell.atoms[jt].tau[j];
				//generates nearest-neighbors shells 
				H_Ewald_pw::rgen(d_tau, rmax, irr, ucell.latvec, ucell.G, r, r2, nrm);
				for(int nr=0 ; nr<nrm ; nr++)
				{
					rr=sqrt(r2[nr]) * ucell.lat0;
					fac = -ModuleBase::e2/2.0/ucell.omega*pow(ucell.lat0,2)*ucell.atoms[it].ncpp.zv * ucell.atoms[jt].ncpp.zv / pow(rr,3) * (erfc(sqa*rr)+rr * sq8a_2pi *  ModuleBase::libm::exp(-alpha * pow(rr,2)));
					for(int l=0; l<3; l++)
					{
						for(int m=0; m<l+1; m++)
						{
							r0[0] = r[nr].x;
							r0[1] = r[nr].y;
							r0[2] = r[nr].z;
							local_sigma(l,m) += fac * r0[l] * r0[m];
						}//end m
					}//end l
				}//end nr
			}

			++ijat;
			ucell.step_jajtiait(&j, &jt, &i, &it);
		}

		delete[] r;
		delete[] r2;
		delete[] irr;
	}//end if

#ifdef _OPENMP
	#pragma omp critical(stress_ewa_reduce)
	{
		sdewald += local_sdewald;
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

	for(int l=0;l<3;l++)
	{
		sigma(l,l) +=sdewald;
	}
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
			sigma(l,m)=-sigma(l,m);
            Parallel_Reduce::reduce_pool(sigma(l, m));
		}
	}
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
			sigma(m,l)=sigma(l,m);
		}
	}

	// this->print(GlobalV::ofs_running, "ewald stress", stression);
	ModuleBase::timer::tick("Stress_Func","stress_ewa");

	return;
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif