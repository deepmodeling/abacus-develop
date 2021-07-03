#include "charge_broyden.h"
#include "global.h"
#include "../src_global/inverse_matrix.h"

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
				for(int is = 0 ; is < NSPIN ; ++is)
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
			for (int is=0; is<NSPIN; is++)
			{
				for (int i=0; i<rstep; i++)
				{
					delete[] Rrho[is][i];
				}
				delete[] Rrho[is];
			}

			// delete: dRrho[i] = Rrho[i+1] - Rrho[i]
			for (int is=0; is<NSPIN; is++)
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

			// dimension of rho_save2(NSPIN, pw.nrxx)
			for (int is=0; is<NSPIN; is++)
			{
				delete[] rho_save2[is];
			}
			delete[] rho_save2;

			// dimension (NSPIN, dstep, dstep)
			delete[] Zmk;

			// dimension (NSPIN, dstep-1, dstep-1)
//		delete[] Zmk_old;
		}
	}
}

void Charge_Broyden::mix_rho
(
    double &dr2,
    const double &diago_error,
    const double &tr2,
    const int &iter,
    bool &converged
)
{
    TITLE("Charge_Broyden","mix_rho");
    timer::tick("Charge","mix_rho",'E');

    for (int is=0; is<NSPIN; is++)
    {
		NOTE("Perform FFT on rho(r) to obtain rho(G).");
        this->set_rhog(rho[is], rhog[is]);

		NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
        this->set_rhog(rho_save[is], rhog_save[is]);


		NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        for (int ig=0; ig<pw.ngmc; ig++)
        {
            this->rhog[is][ig] -= this->rhog_save[is][ig];
        }

    }

	NOTE("Calculate the norm of the Residual vector: < R[rho] | R[rho_save] >");
    dr2 = this->rhog_dot_product( this->rhog, this->rhog);
	
	if(test_charge)ofs_running << " dr2 from rhog_dot_product is " << dr2 << endl;

	// dr2 calculated from real space.
	double dr22 = 0.0;
	for(int is=0; is<NSPIN; is++)
	{
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			dr22 += abs( rho[is][ir] - rho_save[is][ir] );
		}
	}

	Parallel_Reduce::reduce_double_pool( dr22 );
	assert( nelec != 0);
	assert( ucell.omega > 0);
	assert( pw.ncxyz > 0);
	dr22 *= ucell.omega / static_cast<double>( pw.ncxyz );
	dr22 /= nelec;
	if(test_charge)ofs_running << " dr2 from real space grid is " << dr22 << endl;

	// mohan add 2011-01-22
	//if(LINEAR_SCALING && LOCAL_BASIS) xiaohui modify 2013-09-01
	if(BASIS_TYPE=="lcao" )
	{
		dr2 = dr22;	
	}
    if ( dr2 < diago_error )
    {
        ofs_warning << " dr2 < diago_error, keep charge density unchanged." << endl;
    	timer::tick("Charge","mix_rho");
        return;
    }
    else if (dr2 < tr2)
    {
        converged = true;
    	timer::tick("Charge","mix_rho");
		return;
    }


	// the charge before mixing.
	double **rho123 = new double*[NSPIN];
	for(int is=0; is<NSPIN; ++is)
	{
		rho123[is] = new double[pw.nrxx];
		ZEROS(rho123[is], pw.nrxx);
		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			rho123[is][ir] = this->rho[is][ir];
		}
	}

	
	if ( this->mixing_mode == "plain")
    {
        // calculate mixing change, and save it in rho1.
        for (int is=0; is<NSPIN; is++)
        {
            this->plain_mixing( this->rho[is], this->rho_save[is]);
        }
    }
    else if ( this->mixing_mode == "kerker")
    {
        for (int is=0; is<NSPIN; is++)
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
        WARNING_QUIT("Charge_Pulay","Not implemended yet,coming soon.");
    }

	// mohan add 2011-06-07
	//this->renormalize_rho();

	// mohan add 2012-06-05
	// rho_save is the charge before mixing
	for(int is=0; is<NSPIN; is++)
	{
		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			rho_save[is][ir] = rho123[is][ir];
		}
    }

//	for(int is=0; is<NSPIN; ++is)
  //	DCOPY(rho[is],rho_save[is],pw.nrxx);
	//2014-06-22
	for(int is=0; is<NSPIN; ++is)
	{
		delete[] rho123[is];
	}
	delete[] rho123;

    timer::tick("Charge","mix_rho",'E');
    return;
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
		for(int is=0; is<NSPIN; is++)
		{
			for(int ig = 0 ; ig < pw.ngmc; ++ig)
			{
				dF[ipos][is][ig] -= this->rhog[is][ig];
				dn[ipos][is][ig] -= this->rhog_save[is][ig];
			}
		}
	}
	for(int is=0; is<NSPIN; is++)
	{
		for(int ig = 0 ; ig < pw.ngmc; ++ig)
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
		if(info != 0) WARNING_QUIT("Broyden_mixing", "Error when factorizing beta.");
		dsytri_(&uu,&iter_used,beta.c,&iter_used,iwork,work,&info);
		if(info != 0) WARNING_QUIT("Broyden_mixing", "Error when DSYTRI beta.");
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
			for(int is=0; is<NSPIN; is++)
			{
				for(int ig = 0 ; ig < pw.ngmc; ++ig)
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
	
	for(int is=0; is<NSPIN; is++)
	{
		for(int ig = 0 ; ig < pw.ngmc; ++ig)
		{
			dF[inext][is][ig] = dF[mixing_ndim][is][ig];
			dn[inext][is][ig] = dn[mixing_ndim][is][ig];
		}
	}

	//kerker part if needed
	{
	}

	for(int is=0; is<NSPIN; is++)
	{
		for(int ig = 0 ; ig < pw.ngmc; ++ig)
		{
			rhog_save[is][ig] += mixing_beta * rhog[is][ig];
		}
		this->set_rhor( rhog_save[is], rho[is]);
	}



	return;
}

void Charge_Broyden::Modified_Broyden_mixing(void)
{
    //TITLE("Charge_Broyden","Modified_Broyden_Mixing");

    this->rstep = this->mixing_ndim;
    this->dstep = this->rstep - 1;

	//cout << "\n initb = " << initb << endl;
    
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

	//cout << "\n irstep = " << irstep;
	//cout << "\n idstep = " << idstep;
	//cout << "\n totstep = " << totstep;

	this->generate_datas(irstep, idstep, totstep);

	// if not enough step, take kerker mixing method.
	if(totstep < dstep)
	{
		for(int is=0; is<NSPIN; is++)
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
		for(int is=0; is<NSPIN; is++)
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
			this->dF = new complex<double>**[npdim];
    		this->dn = new complex<double>**[npdim];
			
			for (int i=0; i<npdim; i++)
    		{
				dF[i] = new complex<double>*[NSPIN]; 
    	    	dn[i] = new complex<double>*[NSPIN]; 
				for (int is=0; is<NSPIN; is++)
    	    	{
    	        	dF[i][is] = new complex<double>[pw.ngmc];
    	        	dn[i][is] = new complex<double>[pw.ngmc];
    	    	}
			}
			Memory::record("Charge_Broyden","dF", NSPIN*npdim*pw.ngmc,"cdouble");
    		Memory::record("Charge_Broyden","dn", NSPIN*npdim*pw.ngmc,"cdouble");
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
    		this->Rrho = new double**[NSPIN];
    		for (int is=0; is<NSPIN; is++)
    		{
    	    	this->Rrho[is] = new double*[rstep];
    	    	for (int i=0; i<rstep; i++)
    	    	{
    	        	this->Rrho[is][i] = new double[pw.nrxx];
    	        	ZEROS(Rrho[is][i],pw.nrxx);
    	    	}	
    		}
    		Memory::record("Charge_Broyden","Rrho", NSPIN*rstep*pw.nrxx,"double");

    		// (2) allocate dRrho[i]: Rrho[i+1] - Rrho[i]
    		this->dRrho = new double**[NSPIN];
    		this->drho = new double**[NSPIN];
    		this->rho_save2 = new double*[NSPIN];
    		for (int is=0; is<NSPIN; is++)
    		{
    	    	dRrho[is] = new double*[dstep];
    	    	drho[is] = new double*[dstep];
    	    	rho_save2[is] = new double[pw.nrxx];
    	    	for (int i=0; i<dstep; i++)
    	    	{
    	        	dRrho[is][i] = new double[pw.nrxx];
    	        	drho[is][i] = new double[pw.nrxx];
    	        	ZEROS( dRrho[is][i], pw.nrxx );
    	        	ZEROS( drho[is][i], pw.nrxx);
    	    	}
    		}
    		Memory::record("Charge_Broyden","dRrho", NSPIN*dstep*pw.nrxx,"double");
    		Memory::record("Charge_Broyden","drho", NSPIN*dstep*pw.nrxx,"double");
    		Memory::record("Charge_Broyden","rho_save2", NSPIN*pw.nrxx,"double");

			this->dRR = new double[dstep];
			ZEROS(dRR, dstep);

			this->beta.create(dstep, dstep);
			this->betabar.create(dstep, dstep);
			this->Abar.create(dstep, dstep);
			

			this->Zmk = new matrix[NSPIN];
//			this->Zmk_old = new matrix[NSPIN];
			for(int is=0; is<NSPIN; is++)
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
	//TITLE("Charge_Broyden","generate_beta");

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
//	cout << endl;
		
	return;
}

void Charge_Broyden::generate_Zmk(const int &totstep, const int &irstep, const int &idstep, const int &is)
{
	//TITLE("Charge_Bryoden","generate_Zmk");
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

	cout << "\n irstep=" << irstep;
	out.printrm("Zmk",Zmk[is],1.0e-15);
	out.printrm("Zmk_old",Zmk_old[is],1.0e-15);
	cout << endl;
	
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
//	TITLE("Charge_Broyden","generate_new_broyden_rho");
	double mixp = this->mixing_beta;

	// gamma save how much 'u' to mix.
	double* gamma = new double[dstep];
	ZEROS(gamma, dstep);
	for(int i=0; i<dstep; i++)
	{
		for(int k=0; k<dstep; k++)
		{
			gamma[i] += this->Zmk[is](k,i) * this->dRR[k];
		}
		//cout << "\n gamma[" << i << "]=" << gamma[i];
	}
	//cout << endl;
	
	for(int ir=0; ir<pw.nrxx; ir++)
	{
		this->rho[is][ir] = this->rho_save[is][ir] + mixp * this->Rrho[is][m][ir];
	
		for(int i=0; i<dstep; i++)
		{
			this->rho[is][ir] -= gamma[i] * ( this->drho[is][i][ir] + mixp * this->dRrho[is][i][ir] );
		}
	}

	/*
	cout << "\n check E: " << endl;
	
	double* rhot = new double[pw.nrxx];
	ZEROS(rhot, pw.nrxx);

	for(int i=0; i<dstep; i++)
	{
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			rhot[ir] = this->drho[is][i][ir] + gamma[i] * this->dRrho[is][i][ir];
		}
		cout << "\n residual_norm = " << this->calculate_residual_norm(rhot,rhot) << endl;
	}	

	BLOCK_HERE("haha");

	delete[] rhot;
	*/

	DCOPY(rho[is], rho_save[is], pw.nrxx);

	delete[] gamma;
	return;
}
