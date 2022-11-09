#include "charge_broyden.h"
#include "global.h"
#include "../module_base/global_variable.h"
#include "../module_base/inverse_matrix.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

Charge_Broyden::Charge_Broyden() 
{
	initb = false;	
}

Charge_Broyden::~Charge_Broyden()
{
    if (initb)
	{
		for (int i=0; i<mixing_ndim+1; ++i)
		{
			for(int is = 0 ; is < GlobalV::NSPIN ; ++is)
			{
				delete[] dF[i][is];
				delete[] dn[i][is];
			}
			delete[] dF[i];
			delete[] dn[i];
		}
		delete[] dF;
		delete[] dn;
	}
}

double Charge_Broyden::get_drho(double** rho, double** rho_save,
	std::complex<double>** rhog, std::complex<double>** rhog_save, const double nelec)
{
	double scf_thr;
	for (int is=0; is<GlobalV::NSPIN; is++)
    {
		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
        GlobalC::rhopw->real2recip(rho[is], rhog[is]);

		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
        GlobalC::rhopw->real2recip(rho_save[is], rhog_save[is]);


		ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        for (int ig=0; ig<GlobalC::rhopw->npw; ig++)
        {
            rhog[is][ig] -= rhog_save[is][ig];
        }

    }

	ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
    scf_thr = this->rhog_dot_product( rhog, rhog);
	
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from rhog_dot_product is " << scf_thr << std::endl;

	// scf_thr calculated from real space.
	double scf_thr2 = 0.0;
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			scf_thr2 += abs( rho[is][ir] - rho_save[is][ir] );
		}
	}

	Parallel_Reduce::reduce_double_pool( scf_thr2 );
	assert( nelec != 0);
	assert( GlobalC::ucell.omega > 0);
	assert( GlobalC::rhopw->nxyz > 0);
	scf_thr2 *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );
	scf_thr2 /= nelec;
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from real space grid is " << scf_thr2 << std::endl;

	// mohan add 2011-01-22
	//if(LINEAR_SCALING && LOCAL_BASIS) xiaohui modify 2013-09-01
	if(GlobalV::BASIS_TYPE=="lcao" )
	{
		scf_thr = scf_thr2;	
	}
	return scf_thr;
}

void Charge_Broyden::mix_rho
(
    const int &iter,
	double** rho,
	double** rho_save,
	std::complex<double>** rhog,
	std::complex<double>** rhog_save
)
{
    ModuleBase::TITLE("Charge_Broyden","mix_rho");
	ModuleBase::timer::tick("Charge", "mix_rho");

	// the charge before mixing.
	double **rho123 = new double*[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		rho123[is] = new double[GlobalC::rhopw->nrxx];
		ModuleBase::GlobalFunc::ZEROS(rho123[is], GlobalC::rhopw->nrxx);
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ++ir)
		{
			rho123[is][ir] = rho[is][ir];
		}
	}
	
	if ( this->mixing_mode == "plain")
    {
        // calculate mixing change, and save it in rho1.
        for (int is=0; is<GlobalV::NSPIN; is++)
        {
            this->plain_mixing( rho[is], rho_save[is]);
        }
    }
    else if ( this->mixing_mode == "pulay")
    {
        this->Pulay_mixing(rho, rho_save);
    }
    else if ( this->mixing_mode == "broyden")
    {
		this->Simplified_Broyden_mixing(iter, rho, rho_save, rhog, rhog_save);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Pulay","Not implemended yet,coming soon.");
    }

	// mohan add 2012-06-05
	// rho_save is the charge before mixing
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<GlobalC::rhopw->nrxx; ++ir)
		{
			rho_save[is][ir] = rho123[is][ir];
		}
    }

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		delete[] rho123[is];
	}
	delete[] rho123;

	if(new_e_iteration) new_e_iteration = false;

    ModuleBase::timer::tick("Charge","mix_rho");
    return;
}

void Charge_Broyden::Simplified_Broyden_mixing(const int &iter,
	double** rho,
	double** rho_save,
    std::complex<double>** rhog,
    std::complex<double>** rhog_save)
{
	//It is a simplified modified broyden_mixing method.
	//Ref: D.D. Johnson PRB 38, 12807 (1988)
	//Here the weight w0 of the error of the inverse Jacobian is set to 0 and the weight wn of
	//the error of each previous iteration is set to same.

	// (1)
	this->allocate_Broyden();
	
	int iter_used = min(iter-1, mixing_ndim);
	int ipos = iter-2 - int((iter-2)/mixing_ndim) * mixing_ndim;
	if(iter > 1)
	{
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ig = 0 ; ig < GlobalC::rhopw->npw; ++ig)
			{
				dF[ipos][is][ig] -= rhog[is][ig];
				dn[ipos][is][ig] -= rhog_save[is][ig];
			}
		}
	}
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ig = 0 ; ig < GlobalC::rhopw->npw; ++ig)
		{
			dF[mixing_ndim][is][ig] = rhog[is][ig];
			dn[mixing_ndim][is][ig] = rhog_save[is][ig];
		}
	}
	
	if(iter_used > 0)
	{
		this->beta.create(iter_used, iter_used,false);
		for(int i = 0; i < iter_used; ++i)
		{
			for(int j = i; j < iter_used; ++j)
			{
				beta(i,j) = rhog_dot_product( this->dF[i], this->dF[j] );
				if(j != i)
				{
					beta(j,i)=beta(i,j);
				}
			}
		}
		double * work = new double [iter_used];
		int * iwork = new int [iter_used];
		char uu='U';
		int info;
		dsytrf_(&uu,&iter_used,beta.c,&iter_used,iwork,work,&iter_used,&info);
		if(info != 0) ModuleBase::WARNING_QUIT("Broyden_mixing", "Error when factorizing beta.");
		dsytri_(&uu,&iter_used,beta.c,&iter_used,iwork,work,&info);
		if(info != 0) ModuleBase::WARNING_QUIT("Broyden_mixing", "Error when DSYTRI beta.");
		for(int i = 0; i < iter_used; ++i)
		{
			for(int j = i + 1; j < iter_used; ++j)
			{
				beta(i,j) = beta(j,i);
			}
		}
		for(int i = 0 ; i < iter_used ; ++i)
		{
			work[i] = rhog_dot_product( this->dF[i], rhog );
		}
		for(int i = 0 ; i < iter_used ; ++i)
		{
			double gamma0 = 0;
			for(int j = 0; j < iter_used ; ++j)
			{
				gamma0 += beta(i,j) * work[j];
			}
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				for(int ig = 0 ; ig < GlobalC::rhopw->npw; ++ig)
				{
					rhog[is][ig] -= gamma0 * dF[i][is][ig];
					rhog_save[is][ig] -= gamma0 * dn[i][is][ig];
				}
			}
			
		}
		delete[] work;
		delete[] iwork;
	}
	int inext = iter-1 - int((iter-1)/mixing_ndim) * mixing_ndim;
	
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ig = 0 ; ig < GlobalC::rhopw->npw; ++ig)
		{
			dF[inext][is][ig] = dF[mixing_ndim][is][ig];
			dn[inext][is][ig] = dn[mixing_ndim][is][ig];
		}
	}


	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ig = 0 ; ig < GlobalC::rhopw->npw; ++ig)
		{
			rhog_save[is][ig] += mixing_beta * rhog[is][ig];
		}
		GlobalC::rhopw->recip2real( rhog_save[is], rho[is]);
	}

	return;
}

void Charge_Broyden::allocate_Broyden()
{
	if(!initb)
	{
		int npdim = mixing_ndim + 1; // another array is used for temporarily store
		this->dF = new std::complex<double>**[npdim];
		this->dn = new std::complex<double>**[npdim];
		
		for (int i=0; i<npdim; i++)
		{
			dF[i] = new std::complex<double>*[GlobalV::NSPIN]; 
			dn[i] = new std::complex<double>*[GlobalV::NSPIN]; 
			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				dF[i][is] = new std::complex<double>[GlobalC::rhopw->npw];
				dn[i][is] = new std::complex<double>[GlobalC::rhopw->npw];
			}
		}
		ModuleBase::Memory::record("Charge_Broyden","dF", GlobalV::NSPIN*npdim*GlobalC::rhopw->npw,"cdouble");
		ModuleBase::Memory::record("Charge_Broyden","dn", GlobalV::NSPIN*npdim*GlobalC::rhopw->npw,"cdouble");
		this->initb = true;
	}

    return;
}