#include "charge_mixing.h"
#include "global.h"
#include "../module_base/global_variable.h"
#include "../module_base/inverse_matrix.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

void Charge_Mixing::Simplified_Broyden_mixing(const int &iter,
	double** rho,
	double** rho_save,
    std::complex<double>** rhog,
    std::complex<double>** rhog_save)
{
	ModuleBase::TITLE("Charge_Mixing","Simplified_Broyden_mixing");
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
		if(info != 0) ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when factorizing beta.");
		dsytri_(&uu,&iter_used,beta.c,&iter_used,iwork,work,&info);
		if(info != 0) ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when DSYTRI beta.");
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

void Charge_Mixing::allocate_Broyden()
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
		ModuleBase::Memory::record("Charge_Mixing","dF", GlobalV::NSPIN*npdim*GlobalC::rhopw->npw,"cdouble");
		ModuleBase::Memory::record("Charge_Mixing","dn", GlobalV::NSPIN*npdim*GlobalC::rhopw->npw,"cdouble");
		this->initb = true;
	}

    return;
}

void Charge_Mixing::deallocate_Broyden()
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