#include "stress_func.h"
#include "source_hamilt/module_ewald/H_Ewald_pw.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "source_base/libm/libm.h"
#include "source_pw/module_pwdft/global.h"

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
    ModuleBase::TITLE("Stress","stress_ewa");
    ModuleBase::timer::tick("Stress","stress_ewa");

    FPTYPE charge=0;
    for(int it=0; it < ucell.ntype; it++)
	{
		charge = charge + ucell.atoms[it].ncpp.zv * ucell.atoms[it].na;
	}
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    FPTYPE alpha=2.9;
    FPTYPE upperbound=0.0;

	do{
		alpha-=0.1;
		if(alpha==0.0)
		{
			ModuleBase::WARNING_QUIT("stres_ew", "optimal alpha not found");
		}
		upperbound =ModuleBase::e2 * pow(charge,2) * 
         sqrt( 2 * alpha / (ModuleBase::TWO_PI)) 
         * erfc(sqrt(ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
	}
    while(upperbound>1e-7);

    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    FPTYPE sdewald=0.0;
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
	if (PARAM.globalv.gamma_only_pw && is_pw) 
	{
		fact=2.0;
	}
//    else fact=1.0;

#pragma omp parallel
{
	ModuleBase::matrix local_sigma(3, 3);
	FPTYPE local_sdewald = 0;

    FPTYPE g2,g2a;
    FPTYPE arg;
    std::complex<FPTYPE> rhostar;
    FPTYPE sewald;

	#pragma omp for
    for(int ig = 0; ig < rho_basis->npw; ig++)
	{
		if(ig == ig0)  
		{
			continue;
		}
		g2 = rho_basis->gg[ig]* ucell.tpiba2;
		g2a = g2 /4.0/alpha;
		rhostar=std::complex<FPTYPE>(0.0,0.0);

		for(int it=0; it < ucell.ntype; it++)
		{
			for(int i=0; i<ucell.atoms[it].na; i++)
			{
				arg = (rho_basis->gcar[ig] * ucell.atoms[it].tau[i]) * (ModuleBase::TWO_PI);
				FPTYPE sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
				rhostar = rhostar + std::complex<FPTYPE>(ucell.atoms[it].ncpp.zv * cosp,ucell.atoms[it].ncpp.zv * sinp);
			}
		}
		rhostar /= ucell.omega;
		sewald = fact* (ModuleBase::TWO_PI) * ModuleBase::e2 * ModuleBase::libm::exp(-g2a) / g2 * pow(std::abs(rhostar),2);
		local_sdewald -= sewald;
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				local_sigma(l, m) += sewald * ucell.tpiba2 * 2.0 
					* rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / g2 * (g2a + 1);
			}
		}
	}

    //R-space sum here (only for the processor that contains G=0) 
    int mxr = 200;
    int *irr=nullptr;
    ModuleBase::Vector3<FPTYPE> *r;
    FPTYPE *r2=nullptr;
    FPTYPE rr=0.0;
    ModuleBase::Vector3<FPTYPE> d_tau;
    FPTYPE r0[3];
    FPTYPE rmax=0.0;
    int nrm=0;
    FPTYPE fac=0.0;

	if(ig0 >= 0)
	{
		std::vector<ModuleBase::Vector3<FPTYPE>> r(mxr);
		std::vector<FPTYPE> r2(mxr);
		std::vector<int> irr(mxr);

		FPTYPE sqa = sqrt(alpha);
		FPTYPE sq8a_2pi = sqrt(8 * alpha / (ModuleBase::TWO_PI));
		rmax = 4.0/sqa/ucell.lat0;

		#pragma omp for
		for(long long ijat = 0; ijat < ucell.nat * ucell.nat; ijat++)
		{
			int it=0;
			int i=0;
			int jt=0;
			int j=0;
			ucell.ijat2iaitjajt(ijat, &i, &it, &j, &jt);
			if (ucell.atoms[it].na != 0 && ucell.atoms[jt].na != 0)
			{
				//calculate tau[na]-tau[nb]
				d_tau = ucell.atoms[it].tau[i] - ucell.atoms[jt].tau[j];
				//generates nearest-neighbors shells 
				H_Ewald_pw::rgen(d_tau, rmax, irr.data(), ucell.latvec, ucell.G, r.data(), r2.data(), nrm);
				for(int nr=0; nr<nrm; nr++)
				{
					rr=sqrt(r2[nr]) * ucell.lat0;
					fac = -ModuleBase::e2/2.0/ucell.omega*
                          pow(ucell.lat0,2)*ucell.atoms[it].ncpp.zv * ucell.atoms[jt].ncpp.zv 
                          / pow(rr,3) * (erfc(sqa*rr)+rr * sq8a_2pi *  ModuleBase::libm::exp(-alpha * pow(rr,2)));

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
		}
	}//end if

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

	ModuleBase::timer::tick("Stress","stress_ewa");

	return;
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif
