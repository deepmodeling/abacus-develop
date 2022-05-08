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
		if(broyden_type==0)
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
		else
		{
			// delete: Rrho[i] = rho_out[i] - rho_in[i];
			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				for (int i=0; i<rstep; i++)
				{
					delete[] Rrho[is][i];
				}
				delete[] Rrho[is];
			}

			// delete: dRrho[i] = Rrho[i+1] - Rrho[i]
			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				for (int i=0; i<dstep; i++)
				{
					delete[] dRrho[is][i];
					delete[] drho[is][i];
				}
				delete[] dRrho[is];
				delete[] drho[is];
			}
			delete[] dRrho;
			delete[] drho;

			// dimension : dstep
			delete[] w;
			delete[] dRR;

			// dimension of rho_save2(GlobalV::NSPIN, GlobalC::pw.nrxx)
			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				delete[] rho_save2[is];
			}
			delete[] rho_save2;

			// dimension (GlobalV::NSPIN, dstep, dstep)
			delete[] Zmk;

			// dimension (GlobalV::NSPIN, dstep-1, dstep-1)
//		delete[] Zmk_old;
		}
	}
}

double Charge_Broyden::get_drho()
{
	double scf_thr;
	for (int is=0; is<GlobalV::NSPIN; is++)
    {
		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
        this->set_rhog(rho[is], rhog[is]);

		ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
        this->set_rhog(rho_save[is], rhog_save[is]);


		ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
        {
            this->rhog[is][ig] -= this->rhog_save[is][ig];
        }

    }

	ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
    scf_thr = this->rhog_dot_product( this->rhog, this->rhog);
	
	if(GlobalV::test_charge)GlobalV::ofs_running << " scf_thr from rhog_dot_product is " << scf_thr << std::endl;

	// scf_thr calculated from real space.
	double scf_thr2 = 0.0;
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			scf_thr2 += abs( rho[is][ir] - rho_save[is][ir] );
		}
	}

	Parallel_Reduce::reduce_double_pool( scf_thr2 );
	assert( nelec != 0);
	assert( GlobalC::ucell.omega > 0);
	assert( GlobalC::pw.ncxyz > 0);
	scf_thr2 *= GlobalC::ucell.omega / static_cast<double>( GlobalC::pw.ncxyz );
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
    const int &iter
)
{
    ModuleBase::TITLE("Charge_Broyden","mix_rho");
	ModuleBase::timer::tick("Charge", "mix_rho");

	// the charge before mixing.
	double **rho123 = new double*[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		rho123[is] = new double[GlobalC::pw.nrxx];
		ModuleBase::GlobalFunc::ZEROS(rho123[is], GlobalC::pw.nrxx);
		for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
		{
			rho123[is][ir] = this->rho[is][ir];
		}
	}
	
	
	if ( this->mixing_mode == "plain")
    {
        // calculate mixing change, and save it in rho1.
        for (int is=0; is<GlobalV::NSPIN; is++)
        {
            this->plain_mixing( this->rho[is], this->rho_save[is]);
        }
    }
    else if ( this->mixing_mode == "kerker")
    {
        for (int is=0; is<GlobalV::NSPIN; is++)
        {
            this->Kerker_mixing( this->rho[is], this->rhog[is], this->rho_save[is] );
        }
    }
    else if ( this->mixing_mode == "pulay")
    {
        this->Pulay_mixing();
    }
    else if ( this->mixing_mode == "pulay-kerker")//2015-06-15
    {
        this->Pulay_mixing();
    }
    else if ( this->mixing_mode == "broyden")
    {
		this->Simplified_Broyden_mixing(iter);
        //this->Modified_Broyden_mixing();
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Pulay","Not implemended yet,coming soon.");
    }

	// mohan add 2011-06-07
	//this->renormalize_rho();

	// mohan add 2012-06-05
	// rho_save is the charge before mixing
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
		{
			rho_save[is][ir] = rho123[is][ir];
		}
    }

//	for(int is=0; is<GlobalV::NSPIN; ++is)
  //	ModuleBase::GlobalFunc::DCOPY(rho[is],rho_save[is],GlobalC::pw.nrxx);
	//2014-06-22
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		delete[] rho123[is];
	}
	delete[] rho123;

    ModuleBase::timer::tick("Charge","mix_rho");
    return;
}

void Charge_Broyden::tmp_mixrho
(
    double &scf_thr,
    const double &diago_error,
    const double &tr2,
    const int &iter,
    bool &converged
)
{
	scf_thr = get_drho();
	if ( scf_thr < diago_error )
    {
        GlobalV::ofs_warning << " scf_thr < diago_error, keep charge density unchanged." << std::endl;
    	ModuleBase::timer::tick("Charge","mix_rho");
        return;
    }
    else if (scf_thr < tr2)
    {
        converged = true;
    	ModuleBase::timer::tick("Charge","mix_rho");
		return;
    }
	mix_rho(iter);
}

void Charge_Broyden::Simplified_Broyden_mixing(const int &iter)
{
	//It is a simplified modified broyden_mixing method.
	//Ref: D.D. Johnson PRB 38, 12807 (1988)
	//Here the weight w0 of the error of the inverse Jacobian is set to 0 and the weight wn of
	//the error of each previous iteration is set to same.

	// (1)
	this->broyden_type=0;
	this->allocate_Broyden();
	
	int iter_used = min(iter-1, mixing_ndim);
	int ipos = iter-2 - int((iter-2)/mixing_ndim) * mixing_ndim;
	if(iter > 1)
	{
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ig = 0 ; ig < GlobalC::pw.ngmc; ++ig)
			{
				dF[ipos][is][ig] -= this->rhog[is][ig];
				dn[ipos][is][ig] -= this->rhog_save[is][ig];
			}
		}
	}
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ig = 0 ; ig < GlobalC::pw.ngmc; ++ig)
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
			work[i] = rhog_dot_product( this->dF[i], this->rhog );
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
				for(int ig = 0 ; ig < GlobalC::pw.ngmc; ++ig)
				{
					this->rhog[is][ig] -= gamma0 * dF[i][is][ig];
					this->rhog_save[is][ig] -= gamma0 * dn[i][is][ig];
				}
			}
			
		}
		delete[] work;
		delete[] iwork;
	}
	int inext = iter-1 - int((iter-1)/mixing_ndim) * mixing_ndim;
	
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ig = 0 ; ig < GlobalC::pw.ngmc; ++ig)
		{
			dF[inext][is][ig] = dF[mixing_ndim][is][ig];
			dn[inext][is][ig] = dn[mixing_ndim][is][ig];
		}
	}

	//kerker part if needed
	{
	}

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ig = 0 ; ig < GlobalC::pw.ngmc; ++ig)
		{
			rhog_save[is][ig] += mixing_beta * rhog[is][ig];
		}
		this->set_rhor( rhog_save[is], rho[is]);
	}



	return;
}

void Charge_Broyden::Modified_Broyden_mixing(void)
{
    //ModuleBase::TITLE("Charge_Broyden","Modified_Broyden_Mixing");

    this->rstep = this->mixing_ndim;
    this->dstep = this->rstep - 1;

	//std::cout << "\n initb = " << initb << std::endl;
    
	// (1)
	this->broyden_type=1;
	this->allocate_Broyden();
	
	// irstep: iteration step for rstep (Rrho)
	// icstep: iteration step for dstep (dRrho)
	// totstep only used for the first few iterations.
	static int irstep = 0; // count step for rstep
	static int idstep = 0; // coutn step for dstep
	static int totstep = 0;

	if (irstep==rstep) irstep=0;
	if (idstep==dstep) idstep=0;

	//std::cout << "\n irstep = " << irstep;
	//std::cout << "\n idstep = " << idstep;
	//std::cout << "\n totstep = " << totstep;

	this->generate_datas(irstep, idstep, totstep);

	// if not enough step, take kerker mixing method.
	if(totstep < dstep)
	{
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->Kerker_mixing( this->rho[is], this->rhog[is], this->rho_save[is] );
		}
		++irstep;
		++idstep;
		++totstep;
		return;
	}
	else
	
	{
		/*
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			// irstep is m, 
			this->generate_beta(is);
			this->generate_Zmk(totstep, irstep, idstep, is);				
			this->generate_dRR(irstep);
			this->generate_new_broyden_rho(is,irstep);
		}
		*/
		// need to be update like pulay mixing.		
	}
	
	++irstep;
	++idstep;
	++totstep;
    return;
}


void Charge_Broyden::allocate_Broyden()
{
	if(!initb)
	{
		if(broyden_type==0)
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
    	        	dF[i][is] = new std::complex<double>[GlobalC::pw.ngmc];
    	        	dn[i][is] = new std::complex<double>[GlobalC::pw.ngmc];
    	    	}
			}
			ModuleBase::Memory::record("Charge_Broyden","dF", GlobalV::NSPIN*npdim*GlobalC::pw.ngmc,"cdouble");
    		ModuleBase::Memory::record("Charge_Broyden","dn", GlobalV::NSPIN*npdim*GlobalC::pw.ngmc,"cdouble");
		}
		else
		{
    		assert(rstep > 0);

			// weight
    		this->w0 = 0.01;
    		this->w = new double[rstep];
    		for (int i=0; i<rstep; i++)
    		{
    	    	w[i] = 1.0/static_cast<double>(rstep);
    		}

			//special test: Pulay mixing
			//this->w0 = 0.00;
			//for(int i=0; i<rstep; i++)
			//{
			//	w[i] = 1.0;
			//}

			// R[rho_in] = rho_out - rho_in
    		this->Rrho = new double**[GlobalV::NSPIN];
    		for (int is=0; is<GlobalV::NSPIN; is++)
    		{
    	    	this->Rrho[is] = new double*[rstep];
    	    	for (int i=0; i<rstep; i++)
    	    	{
    	        	this->Rrho[is][i] = new double[GlobalC::pw.nrxx];
    	        	ModuleBase::GlobalFunc::ZEROS(Rrho[is][i],GlobalC::pw.nrxx);
    	    	}	
    		}
    		ModuleBase::Memory::record("Charge_Broyden","Rrho", GlobalV::NSPIN*rstep*GlobalC::pw.nrxx,"double");

    		// (2) allocate dRrho[i]: Rrho[i+1] - Rrho[i]
    		this->dRrho = new double**[GlobalV::NSPIN];
    		this->drho = new double**[GlobalV::NSPIN];
    		this->rho_save2 = new double*[GlobalV::NSPIN];
    		for (int is=0; is<GlobalV::NSPIN; is++)
    		{
    	    	dRrho[is] = new double*[dstep];
    	    	drho[is] = new double*[dstep];
    	    	rho_save2[is] = new double[GlobalC::pw.nrxx];
    	    	for (int i=0; i<dstep; i++)
    	    	{
    	        	dRrho[is][i] = new double[GlobalC::pw.nrxx];
    	        	drho[is][i] = new double[GlobalC::pw.nrxx];
    	        	ModuleBase::GlobalFunc::ZEROS( dRrho[is][i], GlobalC::pw.nrxx );
    	        	ModuleBase::GlobalFunc::ZEROS( drho[is][i], GlobalC::pw.nrxx);
    	    	}
    		}
    		ModuleBase::Memory::record("Charge_Broyden","dRrho", GlobalV::NSPIN*dstep*GlobalC::pw.nrxx,"double");
    		ModuleBase::Memory::record("Charge_Broyden","drho", GlobalV::NSPIN*dstep*GlobalC::pw.nrxx,"double");
    		ModuleBase::Memory::record("Charge_Broyden","rho_save2", GlobalV::NSPIN*GlobalC::pw.nrxx,"double");

			this->dRR = new double[dstep];
			ModuleBase::GlobalFunc::ZEROS(dRR, dstep);

			this->beta.create(dstep, dstep);
			this->betabar.create(dstep, dstep);
			this->Abar.create(dstep, dstep);
			

			this->Zmk = new ModuleBase::matrix[GlobalV::NSPIN];
//			this->Zmk_old = new matrix[GlobalV::NSPIN];
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				this->Zmk[is].create(dstep, dstep);
//				this->Zmk_old[is].create(dstep, dstep);
			}
		}
		this->initb = true;
	}

    return;
}

void Charge_Broyden::generate_beta(const int &is)
{
	//ModuleBase::TITLE("Charge_Broyden","generate_beta");

	//(1) generate Abar(k,n) = w(k)*w(n)*<dR(n)|dR(k)>
	for(int k=0; k<dstep; k++)
	{
		for(int n=0; n<dstep; n++)
		{
			this->Abar(k,n) = this->w[k]*this->w[n]*
			this->calculate_residual_norm( this->dRrho[is][n], this->dRrho[is][k] );
			this->Abar(n,k) = this->Abar(k,n);
		}
	}

	//out.printrm("Abar",Abar,1.0e-15);

	//(2) generate beta(k,n)=(w0*w0+Abar)(k,n)^-1
	for(int k=0; k<dstep; k++)
	{
		for(int n=0; n<dstep; n++)
		{
			this->beta(k,n) = this->Abar(k,n);
			if(n==k) this->beta(k,n) += this->w0 * this->w0;
		}
	}

	// scheme == 1: 
	const int scheme = 1;
	this->inverse_real_symmetry_matrix(scheme,this->beta);
	
	//(3) generate betabar
	for(int k=0; k<dstep; k++)
	{
		for(int n=0; n<dstep; n++)
		{
			if(k==n) this->betabar(k,n) = 1.0;
			else this->betabar(k,n) = 0.0;
			for(int j=0; j<dstep; j++)
			{
				this->betabar(k,n) -= w[k]*w[j]*beta(k,j)*(Abar(n,j)/w[n]*w[j]);
			}
		}
	}
	
//	out.printrm("beta",beta,1.0e-15);
//	out.printrm("betabar",betabar,1.0e-15);
//	std::cout << std::endl;
		
	return;
}

void Charge_Broyden::generate_Zmk(const int &totstep, const int &irstep, const int &idstep, const int &is)
{
	//ModuleBase::TITLE("Charge_Bryoden","generate_Zmk");
	this->Zmk[is].zero_out();
		
	for(int k=0; k<dstep; k++)
	{
		for(int n=0; n<dstep; n++)
		{
			this->Zmk[is](k,n) = this->beta(k,n)*w[k]*w[n];;
		}
	}
//	out.printrm("Zmk",Zmk[is],1.0e-15);
//	out.printrm("Zmk_old",Zmk_old[is],1.0e-15);

	/*	
	for(int k=0; k<dstep; k++)
	{		
		// Zmk = sum( betabar(k,n) * Zmk_old(n) )
 		for(int n=0; n<dstep; n++)
		{
			// only do (dstep-1) step.
			if(n==irstep)continue;
			for(int nn=0; nn<dstep; nn++)
			{
				if(nn==irstep)continue;
				this->Zmk[is](k,n) += this->betabar(k,nn) * this->Zmk_old[is](nn,n);
			}
		}
	}

	std::cout << "\n irstep=" << irstep;
	out.printrm("Zmk",Zmk[is],1.0e-15);
	out.printrm("Zmk_old",Zmk_old[is],1.0e-15);
	std::cout << std::endl;
	
	// save Zmk old
	// sacrifice liite memory to make coding convenient!
	for(int i=0; i<dstep; i++)
	{
		for(int j=0; j<dstep; j++)
		{
			this->Zmk_old[is](i,j) = this->Zmk[is](i,j);
		}
	}

	*/	

	return;
}


void Charge_Broyden::generate_new_broyden_rho(const int &is, const int &m)
{
//	ModuleBase::TITLE("Charge_Broyden","generate_new_broyden_rho");
	double mixp = this->mixing_beta;

	// gamma save how much 'u' to mix.
	double* gamma = new double[dstep];
	ModuleBase::GlobalFunc::ZEROS(gamma, dstep);
	for(int i=0; i<dstep; i++)
	{
		for(int k=0; k<dstep; k++)
		{
			gamma[i] += this->Zmk[is](k,i) * this->dRR[k];
		}
		//std::cout << "\n gamma[" << i << "]=" << gamma[i];
	}
	//std::cout << std::endl;
	
	for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		this->rho[is][ir] = this->rho_save[is][ir] + mixp * this->Rrho[is][m][ir];
	
		for(int i=0; i<dstep; i++)
		{
			this->rho[is][ir] -= gamma[i] * ( this->drho[is][i][ir] + mixp * this->dRrho[is][i][ir] );
		}
	}

	/*
	std::cout << "\n check E: " << std::endl;
	
	double* rhot = new double[GlobalC::pw.nrxx];
	ModuleBase::GlobalFunc::ZEROS(rhot, GlobalC::pw.nrxx);

	for(int i=0; i<dstep; i++)
	{
		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			rhot[ir] = this->drho[is][i][ir] + gamma[i] * this->dRrho[is][i][ir];
		}
		std::cout << "\n residual_norm = " << this->calculate_residual_norm(rhot,rhot) << std::endl;
	}	

	BLOCK_HERE("haha");

	delete[] rhot;
	*/

	ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is], GlobalC::pw.nrxx);

	delete[] gamma;
	return;
}
