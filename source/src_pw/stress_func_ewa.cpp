#include "./stress_func.h"
#include "./H_Ewald_pw.h"
#include "../module_base/timer.h"
#include "global.h"

//calcualte the Ewald stress term in PW and LCAO
void Stress_Func::stress_ewa(ModuleBase::matrix& sigma, ModulePW::PW_Basis* rho_basis, const bool is_pw)
{
    ModuleBase::TITLE("Stress_Func","stress_ewa");
    ModuleBase::timer::tick("Stress_Func","stress_ewa");

    double charge=0;
    for(int it=0; it < GlobalC::ucell.ntype; it++)
	{
		for(int i=0; i<GlobalC::ucell.atoms[it].na; i++)
		{
			charge = charge + GlobalC::ucell.atoms[it].ncpp.zv;
		}
	}
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    double alpha=2.9;
    double upperbound;
    do{
       alpha-=0.1;
       if(alpha==0.0)
          ModuleBase::WARNING_QUIT("stres_ew", "optimal alpha not found");
       upperbound =ModuleBase::e2 * pow(charge,2) * sqrt( 2 * alpha / (ModuleBase::TWO_PI)) * erfc(sqrt(GlobalC::ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    }
    while(upperbound>1e-7);

    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    double sdewald;
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

    double fact=1.0;
    if (INPUT.gamma_only && is_pw) fact=2.0;
//    else fact=1.0;

    double g2,g2a;
    double arg;
    std::complex<double> rhostar;
    double sewald;
    for(int ig = 0; ig < rho_basis->npw; ig++)
	{
		if(ig == ig0)  continue;
		g2 = rho_basis->gg[ig]* GlobalC::ucell.tpiba2;
		g2a = g2 /4.0/alpha;
		rhostar=std::complex<double>(0.0,0.0);
		for(int it=0; it < GlobalC::ucell.ntype; it++)
		{
			for(int i=0; i<GlobalC::ucell.atoms[it].na; i++)
			{
				arg = (rho_basis->gcar[ig] * GlobalC::ucell.atoms[it].tau[i]) * (ModuleBase::TWO_PI);
				rhostar = rhostar + std::complex<double>(GlobalC::ucell.atoms[it].ncpp.zv * cos(arg),GlobalC::ucell.atoms[it].ncpp.zv * sin(arg));
			}
		}
		rhostar /= GlobalC::ucell.omega;
		sewald = fact* (ModuleBase::TWO_PI) * ModuleBase::e2 * exp(-g2a) / g2 * pow(abs(rhostar),2);
		sdewald = sdewald - sewald;
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				sigma(l, m) += sewald * GlobalC::ucell.tpiba2 * 2.0 * rho_basis->gcar[ig][l] * rho_basis->gcar[ig][m] / g2 * (g2a + 1);
			}
		}
	}

	for(int l=0;l<3;l++)
	{
		sigma(l,l) +=sdewald;
	}
    //R-space sum here (only for the processor that contains G=0) 
    int mxr = 50;
    int *irr;
    ModuleBase::Vector3<double> *r;
    double *r2;
    r  = new ModuleBase::Vector3<double>[mxr];
    r2 = new double[mxr];
    irr = new int[mxr];
    double rr;
    ModuleBase::Vector3<double> d_tau;
    double r0[3];
    double rmax=0.0;
    int nrm=0;
    double fac;
	if(ig0 >= 0)
	{
		rmax = 4.0/sqrt(alpha)/GlobalC::ucell.lat0;
		//with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2*10^-1)
		for(int it=0; it < GlobalC::ucell.ntype; it++)
		{
			for(int i=0; i<GlobalC::ucell.atoms[it].na; i++)
			{
				for(int jt=0; jt < GlobalC::ucell.ntype; jt++)
				{
					for(int j=0; j<GlobalC::ucell.atoms[jt].na; j++)
					{
						//calculate tau[na]-tau[nb]
						d_tau = GlobalC::ucell.atoms[it].tau[i] - GlobalC::ucell.atoms[jt].tau[j];
						//generates nearest-neighbors shells 
						H_Ewald_pw::rgen(d_tau, rmax, irr, GlobalC::ucell.latvec, GlobalC::ucell.G, r, r2, nrm);
						for(int nr=0 ; nr<nrm ; nr++)
						{
							rr=sqrt(r2[nr]) * GlobalC::ucell.lat0;
							fac = -ModuleBase::e2/2.0/GlobalC::ucell.omega*pow(GlobalC::ucell.lat0,2)*GlobalC::ucell.atoms[it].ncpp.zv * GlobalC::ucell.atoms[jt].ncpp.zv / pow(rr,3) * (erfc(sqrt(alpha)*rr)+rr * sqrt(8 * alpha / (ModuleBase::TWO_PI)) * exp(-alpha * pow(rr,2)));
							for(int l=0; l<3; l++)
							{
								for(int m=0; m<l+1; m++)
								{
									r0[0] = r[nr].x;
									r0[1] = r[nr].y;
									r0[2] = r[nr].z;
									sigma(l,m) += fac * r0[l] * r0[m];
								}//end m
							}//end l
						}//end nr
					}//end j
				}//end jt
			}//end i
		}//end it
	}//end if
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
			sigma(m,l)=sigma(l,m);
		}
	}
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigma(l,m)=-sigma(l,m);
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}

	delete[] r;
	delete[] r2;
	delete[] irr;
	// this->print(GlobalV::ofs_running, "ewald stress", stression);
	ModuleBase::timer::tick("Stress_Func","stress_ewa");

	return;
}
