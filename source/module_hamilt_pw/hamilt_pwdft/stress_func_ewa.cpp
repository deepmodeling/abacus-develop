#include "stress_func.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/libm/libm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//calcualte the Ewald stress term in PW and LCAO
template<typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_ewa(ModuleBase::matrix& sigma, ModulePW::PW_Basis* rho_basis, const bool is_pw)
{
    ModuleBase::TITLE("Stress_Func","stress_ewa");
    ModuleBase::timer::tick("Stress_Func","stress_ewa");

    FPTYPE charge=0;
    for(int it=0; it < GlobalC::ucell.ntype; it++)
	{
		charge = charge + GlobalC::ucell.atoms[it].ncpp.zv * GlobalC::ucell.atoms[it].na;
	}
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    FPTYPE alpha=2.9;
    FPTYPE upperbound;
    do{
       alpha-=0.1;
       if(alpha==0.0)
          ModuleBase::WARNING_QUIT("stres_ew", "optimal alpha not found");
       upperbound =ModuleBase::e2 * pow(charge,2) * sqrt( 2 * alpha / (ModuleBase::TWO_PI)) * erfc(sqrt(GlobalC::ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    }
    while(upperbound>1e-7);

    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    FPTYPE sdewald;
	const int ig0 = rho_basis->ig_gge0;
    if( ig0 >= 0)
	{
       sdewald = (ModuleBase::TWO_PI) * ModuleBase::e2 / 4.0 / alpha * pow(charge/GlobalC::ucell.omega,2);
    }
    else 
	{
       sdewald = 0.0;
    }

    //sdewald is the diagonal term 

    FPTYPE fact=1.0;
    if (INPUT.gamma_only && is_pw) fact=2.0;
//    else fact=1.0;

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
    FPTYPE arg;
    std::complex<FPTYPE> rhostar;
    FPTYPE sewald;
    for(; ig < ig_end; ig++)
	{
		if(ig == ig0)  continue;
		g2 = rho_basis->gg[ig]* GlobalC::ucell.tpiba2;
		g2a = g2 /4.0/alpha;
		rhostar=std::complex<FPTYPE>(0.0,0.0);
		for(int it=0; it < GlobalC::ucell.ntype; it++)
		{
			for(int i=0; i<GlobalC::ucell.atoms[it].na; i++)
			{
				arg = (rho_basis->gcar[ig] * GlobalC::ucell.atoms[it].tau[i]) * (ModuleBase::TWO_PI);
				FPTYPE sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
				rhostar = rhostar + std::complex<FPTYPE>(GlobalC::ucell.atoms[it].ncpp.zv * cosp,GlobalC::ucell.atoms[it].ncpp.zv * sinp);
			}
		}
		rhostar /= GlobalC::ucell.omega;
		sewald = fact* (ModuleBase::TWO_PI) * ModuleBase::e2 * ModuleBase::libm::exp(-g2a) / g2 * pow(std::abs(rhostar),2);
		local_sdewald -= sewald;
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				local_sigma(l, m) += sewald * GlobalC::ucell.tpiba2 * 2.0 * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / g2 * (g2a + 1);
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
		rmax = 4.0/sqa/GlobalC::ucell.lat0;

		// collapse it, ia, jt, ja loop into a single loop
		long long ijat, ijat_end;
		int it, i, jt, j;
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, (long long)GlobalC::ucell.nat * GlobalC::ucell.nat, ijat, ijat_end);
		ijat_end = ijat + ijat_end;
		GlobalC::ucell.ijat2iaitjajt(ijat, &i, &it, &j, &jt);

		while (ijat < ijat_end)
		{
			//calculate tau[na]-tau[nb]
			d_tau = GlobalC::ucell.atoms[it].tau[i] - GlobalC::ucell.atoms[jt].tau[j];
			//generates nearest-neighbors shells 
			H_Ewald_pw::rgen(d_tau, rmax, irr, GlobalC::ucell.latvec, GlobalC::ucell.G, r, r2, nrm);
			for(int nr=0 ; nr<nrm ; nr++)
			{
				rr=sqrt(r2[nr]) * GlobalC::ucell.lat0;
				fac = -ModuleBase::e2/2.0/GlobalC::ucell.omega*pow(GlobalC::ucell.lat0,2)*GlobalC::ucell.atoms[it].ncpp.zv * GlobalC::ucell.atoms[jt].ncpp.zv / pow(rr,3) * (erfc(sqa*rr)+rr * sq8a_2pi *  ModuleBase::libm::exp(-alpha * pow(rr,2)));
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

			++ijat;
			GlobalC::ucell.step_jajtiait(&j, &jt, &i, &it);
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

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif