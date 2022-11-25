#include "charge_mixing.h"
#include "global.h"
#include "../module_base/inverse_matrix.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/memory.h"

void Charge_Mixing::Pulay_mixing(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing","Pulay_mixing");
	rstep = this->mixing_ndim;
	dstep = this->mixing_ndim - 1;
	assert(dstep>0);

	// (1) allocate
	this->allocate_Pulay();

    if (irstep==rstep) irstep=0;
	if (idstep==dstep) idstep=0;

	if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"irstep",irstep);
	if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"idstep",idstep);
	if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"totstep",totstep);

	//----------------------------------------------
	// calculate "dR^{i} = R^{i+1} - R^{i}"
	// calculate "drho^{i} = rho^{i+1} - rho^{i}"
	//----------------------------------------------
	this->generate_datas(irstep, idstep, totstep, chr->rho, chr->rho_save);

	// not enough steps, not full matrix.
	if(totstep < dstep) 
	{
		int premix = 1;
		
		// mohan update 2011-06-14
		if(totstep==0) premix = 2;

		if(premix == 1)
		{
			this->generate_Abar(Abar);
			
			//-----------------------------
			// inverse part of the matrix.
			//-----------------------------
			const int predim = idstep+1;
			ModuleBase::matrix preA(predim, predim);
			for(int i=0; i<predim; i++)
			{
				for(int j=0; j<predim; j++)
				{
					preA(i,j) = Abar(i,j);
				}
			}
			this->inverse_preA(predim,preA);
			
			
			//----------------------------------------------
			// get the information from part of the matrix. 
			//----------------------------------------------
			Abar.zero_out();
			for(int i=0; i<predim; i++)
			{
				for(int j=0; j<predim; j++)
				{
					Abar(i,j) = preA(i,j);
				}			
			}

			this->generate_dRR(irstep);

			this->generate_alpha();

			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				this->generate_new_rho(is,irstep,chr->rho,chr->rho_save);
			}
		}

		if(premix == 2)
		{
			// if not enough step, take kerker mixing method.	
			this->plain_mixing(chr);
		}

		if(totstep>0)
		{
			++idstep;
		}
		++irstep;
		++totstep;
	}
	else
	{
		// generate matrix A = <dR|dR>
		this->generate_Abar(Abar);
		ModuleBase::matrix A(Abar);

		// inverse A matrix to become Abar.
		this->inverse_real_symmetry_matrix(Abar);
		
		this->generate_dRR(irstep);

		this->generate_alpha();

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->generate_new_rho(is,irstep,chr->rho,chr->rho_save);
		}

		++irstep;
		++idstep;
		++totstep;	
	}

	return;		
}

void Charge_Mixing::reset()		// Peize Lin add 2018-11-01
{
	this->new_e_iteration = true;
	
	irstep = 0;
	idstep = 0;
	totstep = 0;
}

void Charge_Mixing::allocate_Pulay()
{
	if(!this->initp)
	{
		ModuleBase::TITLE("Charge_Mixing","allocate_pulay");
		ModuleBase::GlobalFunc::NOTE("rstep is used to record Rrho");
		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"rstep",rstep);
		ModuleBase::GlobalFunc::NOTE("dstep is used to record dRrho, drho");
		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dstep",dstep);
		assert(rstep>1);
		dstep = rstep - 1;
		
		// (1) allocate Rrho[i]: rho_out[i] - rho_in[i]
    	this->Rrho = new double**[GlobalV::NSPIN];
    	for (int is=0; is<GlobalV::NSPIN; is++)
    	{
        	Rrho[is] = new double*[rstep];
        	for (int i=0; i<rstep; i++)
        	{
            	Rrho[is][i] = new double[GlobalC::rhopw->nrxx];
				ModuleBase::GlobalFunc::ZEROS( Rrho[is][i], GlobalC::rhopw->nrxx );
        	}
		}
    	ModuleBase::Memory::record("Charge_Mixing","Rrho", GlobalV::NSPIN*rstep*GlobalC::rhopw->nrxx,"double");

		// (2) allocate "dRrho[i] = Rrho[i+1] - Rrho[i]" of the last few steps.
		// allocate "drho[i] = rho[i+1] - rho[i]" of the last few steps.
		// allocate rho_save2: rho[i] of the last few steps.
		this->dRrho = new double**[GlobalV::NSPIN];
		this->drho = new double**[GlobalV::NSPIN];
		this->rho_save2 = new double*[GlobalV::NSPIN];

		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			dRrho[is] = new double*[dstep];
			drho[is] = new double*[dstep];
			rho_save2[is] = new double[GlobalC::rhopw->nrxx];
			ModuleBase::GlobalFunc::ZEROS( rho_save2[is], GlobalC::rhopw->nrxx);

			for (int i=0; i<dstep; i++)
			{
				dRrho[is][i] = new double[GlobalC::rhopw->nrxx];	
				drho[is][i] = new double[GlobalC::rhopw->nrxx];
				ModuleBase::GlobalFunc::ZEROS( dRrho[is][i], GlobalC::rhopw->nrxx );
				ModuleBase::GlobalFunc::ZEROS( drho[is][i], GlobalC::rhopw->nrxx);
			}
		}
    	ModuleBase::Memory::record("Charge_Mixing","dRrho", GlobalV::NSPIN*dstep*GlobalC::rhopw->nrxx,"double");
    	ModuleBase::Memory::record("Charge_Mixing","drho", GlobalV::NSPIN*dstep*GlobalC::rhopw->nrxx,"double");
    	ModuleBase::Memory::record("Charge_Mixing","rho_save2", GlobalV::NSPIN*GlobalC::rhopw->nrxx,"double");

		ModuleBase::GlobalFunc::NOTE("Allocate Abar = <dRrho_j | dRrho_i >, dimension = dstep.");
		this->Abar.create(dstep, dstep);
		ModuleBase::Memory::record("Charge_Mixing","Abar", dstep*dstep,"double");

		// (4) allocate dRR = <delta R|R>
		ModuleBase::GlobalFunc::NOTE("Allocate dRR = < dR | R >, dimension = dstep");
		this->dRR = new double[dstep];
		ModuleBase::GlobalFunc::ZEROS(dRR, dstep);

		// (5) allocate alpha
		ModuleBase::GlobalFunc::NOTE("Allocate alpha, dimension = dstep");
		this->alpha = new double[dstep];
		ModuleBase::GlobalFunc::ZEROS(alpha, dstep);

		this->initp = true;
    }

	// mohan add 2010-07-16
	if(this->new_e_iteration)
	{
		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			for(int j=0; j<rstep; j++)
			{
				ModuleBase::GlobalFunc::ZEROS(Rrho[i][j], GlobalC::rhopw->nrxx);
			}
		}
		
		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			ModuleBase::GlobalFunc::ZEROS(rho_save2[i], GlobalC::rhopw->nrxx);
			for(int j=0; j<dstep; j++)
			{
				ModuleBase::GlobalFunc::ZEROS( dRrho[i][j], GlobalC::rhopw->nrxx );
				ModuleBase::GlobalFunc::ZEROS( drho[i][j], GlobalC::rhopw->nrxx);
			}
		}

		ModuleBase::GlobalFunc::ZEROS(dRR, dstep);
		ModuleBase::GlobalFunc::ZEROS(alpha, dstep);
	}
	return;
}

void Charge_Mixing::deallocate_Pulay()
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
	delete[] Rrho;

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

	// dimension: dstep
	delete[] dRR;
	delete[] alpha;

	// dimension of rho_save2(GlobalV::NSPIN, GlobalC::rhopw->nrxx)
	for (int is=0; is<GlobalV::NSPIN; is++)
	{
		delete[] rho_save2[is];
	}
	delete[] rho_save2;	
}

// calculate < dR | dR >
// if spin is considered, double the size.
// < dR1,dR2 | dR1,dR2 > = < dR1 | dR1 > + < dR2 | dR2 >
void Charge_Mixing::generate_Abar(ModuleBase::matrix &A)const
{
	int step = 0;

	step=dstep;

	A.zero_out();

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int i=0; i<step; i++)
		{
			for(int j=0; j<=i; j++)
			{
				A(i,j) += this->calculate_residual_norm( this->dRrho[is][j], this->dRrho[is][i] );
				A(j,i) = A(i,j);
			}
		}
	}

	return;
}

#include "../module_base/complexmatrix.h"
void Charge_Mixing::inverse_preA(const int &dim, ModuleBase::matrix &preA)const
{
	ModuleBase::ComplexMatrix B(dim, dim);
	ModuleBase::ComplexMatrix C(dim, dim);
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			B(i,j) = std::complex<double> (preA(i,j), 0.0);
		}
	}
	ModuleBase::Inverse_Matrix_Complex IMC;
	IMC.init(dim);
	IMC.using_zheev(B,C);

	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			preA(i,j) = C(i,j).real();
		}
	}
	return;
}

void Charge_Mixing::inverse_real_symmetry_matrix(ModuleBase::matrix &A)const // indicate the spin.
{
//	ModuleBase::TITLE("Charge_Mixing","inverse_Abar");

	int step = 0;

	step=dstep;

	// Notice that it's a symmetry matrix!!!	
	ModuleBase::ComplexMatrix B(step,step);
	ModuleBase::ComplexMatrix C(step,step);
	for(int i=0; i<step; i++)
	{
		for(int j=0; j<step; j++)
		{
			B(i,j) = std::complex<double> (A(i,j),0.0);
		}
	}
		
	ModuleBase::Inverse_Matrix_Complex IMC;
	IMC.init(step);
	IMC.using_zheev(B,C);

	for(int i=0; i<step; i++)
	{
		for(int j=0; j<step; j++)
		{
			A(i,j) = C(i,j).real();
		}
	}
	
	return;
}

void Charge_Mixing::generate_dRR(const int &m)
{
	ModuleBase::GlobalFunc::ZEROS(dRR, dstep);
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int i=0; i<dstep; i++)
		{
			this->dRR[i] += this->calculate_residual_norm(this->dRrho[is][i], this->Rrho[is][m]);	
		}
	}

	return;
}

// use dstep to genearte Abar(dstep, dstep)
void Charge_Mixing::generate_alpha()
{
//	ModuleBase::TITLE("Charge_Mixing","generate_alpha");

	ModuleBase::GlobalFunc::ZEROS(alpha, dstep);
	for(int i=0; i<dstep; i++)
	{
		for(int j=0; j<dstep; j++)
		{
			// Abar is a symmetry matrix.
			this->alpha[i] -= this->Abar(j,i) * this->dRR[j];
		}	
	}

	return;
}

void Charge_Mixing::generate_new_rho(const int &is, const int &m, double** rho, double** rho_save)
{
//	ModuleBase::TITLE("Charge_Mixing","generate_new_rho");
	
	double mixp = this->mixing_beta;
	
	// rho tmp
	double* rhonew = new double[GlobalC::rhopw->nrxx];
	ModuleBase::GlobalFunc::ZEROS(rhonew, GlobalC::rhopw->nrxx);
	
	for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	{
		rhonew[ir] = rho_save[is][ir] + mixp * this->Rrho[is][m][ir];
		for(int i=0; i<dstep; i++)
		{
			rhonew[ir] += this->alpha[i] * ( this->drho[is][i][ir] + mixp * this->dRrho[is][i][ir] );
		}
	}

	ModuleBase::GlobalFunc::DCOPY(rhonew, rho[is], GlobalC::rhopw->nrxx);
	delete[] rhonew;

	return;
}

void Charge_Mixing::generate_residual_vector(double *residual, const double* rho_out, const double* rho_in)const
{
	for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	{
		residual[ir]= rho_out[ir] - rho_in[ir];
	}
	return;
}

double Charge_Mixing::calculate_residual_norm(double *residual1, double* residual2)const
{
	// calculate the norm of the residual std::vector:
	// (the target to minimize in Pulay's algorithm)
	double rnorm = 0.0;
	for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	{
		rnorm += residual1[ir]*residual2[ir];
	}
	Parallel_Reduce::reduce_double_pool(rnorm); // mohan fix bug 2010-07-22
	return rnorm;
}


// calculate "dR^{i} = R^{i+1} - R^{i}"
// calculate "drho^{i} = rho^{i+1} - rho^{i}"
void Charge_Mixing::generate_datas(const int &irstep, const int &idstep, const int &totstep, double** rho, double** rho_save)
{
	//===============================================
	// calculate the important "Rrho".
	// init the "Rrho = rho - rho_save"
	// which Rrho to be update now? answer: irstep
	//===============================================

	ModuleBase::GlobalFunc::NOTE("Generate Residual std::vector from rho and rho_save.");
    for (int is=0; is<GlobalV::NSPIN; is++)
    {
		this->generate_residual_vector( this->Rrho[is][irstep], rho[is], rho_save[is]);

		if(this->mixing_gg0 > 0.0)
		{
			std::complex<double> *kerpulay = new std::complex<double>[GlobalC::rhopw->npw];
			double* kerpulayR = new double[GlobalC::rhopw->nrxx];
			
			GlobalC::rhopw->real2recip(Rrho[is][irstep], kerpulay);

			const double fac = this->mixing_gg0;
			const double gg0 = std::pow(fac * 0.529177 /GlobalC::ucell.tpiba, 2);
			double* filter_g = new double[GlobalC::rhopw->npw];
			for(int ig=0; ig<GlobalC::rhopw->npw; ig++)
			{
				double gg = GlobalC::rhopw->gg[ig];
				filter_g[ig] = max(gg / (gg + gg0), 0.1);
				kerpulay[ig] = (1 - filter_g[ig]) * kerpulay[ig];
			}

			GlobalC::rhopw->recip2real(kerpulay, kerpulayR);
			for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
			{
				Rrho[is][irstep][ir] = Rrho[is][irstep][ir] - kerpulayR[ir];
			}

			delete[] kerpulay;
			delete[] kerpulayR;
			delete[] filter_g;
		}
    }

	if(totstep==0)
	{
		// don't need to calculate 'dRrho' and 'drho'.
	}
	else if(totstep>0)
	{	
		// which dRrho to be update now? answer: dstep=istep-1;
		// irstep range: [0, rstep)
		const int nowR = irstep;
		int lastR = irstep - 1;

		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "now irstep", nowR);
		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "last irstep",lastR); 

		if(lastR < 0) lastR += rstep;
		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
			{
				this->dRrho[is][idstep][ir] = this->Rrho[is][nowR][ir] - this->Rrho[is][lastR][ir];
				this->drho[is][idstep][ir] = rho_save[is][ir] - this->rho_save2[is][ir];
			}
		}
	}

	// mohan fix the bug 2010/03/26.
	// save 'rho_save2' in order to calculate drho in
	// the next iteration.
	ModuleBase::GlobalFunc::NOTE("Calculate drho = rho_{in}^{i+1} - rho_{in}^{i}");
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		ModuleBase::GlobalFunc::DCOPY(rho_save[is], this->rho_save2[is], GlobalC::rhopw->nrxx);
	}
	return;
}

