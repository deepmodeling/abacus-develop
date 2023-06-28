#include "charge_mixing.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/inverse_matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/memory.h"
#include "module_base/tool_threading.h"
#include "module_base/timer.h"

static inline double calculate_residual_norm(int nrxx, double *residual1, double* residual2)
{
	// calculate the norm of the residual std::vector:
	// (the target to minimize in Pulay's algorithm)
	double rnorm = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:rnorm)
#endif
	for(int ir=0; ir<nrxx; ir++)
	{
		rnorm += residual1[ir]*residual2[ir];
	}
	return rnorm;
}

void Charge_Mixing::Pulay_mixing(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing","Pulay_mixing");
	ModuleBase::timer::tick("Charge", "Pulay_mixing");
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
	this->generate_datas(irstep, idstep, totstep, chr);

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
				this->generate_new_rho(is,irstep,chr);
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

		// inverse A matrix to become Abar.
		this->inverse_real_symmetry_matrix(Abar);
		
		this->generate_dRR(irstep);

		this->generate_alpha();

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->generate_new_rho(is,irstep,chr);
		}

		++irstep;
		++idstep;
		++totstep;	
	}

	ModuleBase::timer::tick("Charge", "Pulay_mixing");
	return;		
}

void Charge_Mixing::reset()		// Peize Lin add 2018-11-01
{
	this->new_e_iteration = true;
	
	irstep = 0;
	idstep = 0;
	totstep = 0;

    // liuyu add 2023-03-29
    // if md_prec_level == 2, charge mixing should re-allocate 
    // due to the change of FFT grids
    if (GlobalV::md_prec_level == 2)
    {
        if (this->mixing_mode == "pulay")
        {
            this->deallocate_Pulay();
        }
        else if (this->mixing_mode == "broyden")
        {
            this->deallocate_Broyden();
        }
    }
}

void Charge_Mixing::set_rhopw(ModulePW::PW_Basis* rhopw_in)
{
	this->rhopw = rhopw_in;
}

void Charge_Mixing::allocate_Pulay()
{
	auto zeros_kernel = [&](int num_threads, int thread_id)
	{
		int beg, len;
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, this->rhopw->nrxx, beg, len);
		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			for(int j=0; j<rstep; j++)
			{
				ModuleBase::GlobalFunc::ZEROS(Rrho[i][j] + beg, len);
			}
		}

		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			ModuleBase::GlobalFunc::ZEROS(rho_save2[i] + beg, len);
			for(int j=0; j<dstep; j++)
			{
				ModuleBase::GlobalFunc::ZEROS( dRrho[i][j] + beg, len );
				ModuleBase::GlobalFunc::ZEROS( drho[i][j] + beg, len);
			}
		}

		if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
		{
			for(int i=0; i<GlobalV::NSPIN; i++)
			{
				for(int j=0; j<rstep; j++)
				{
					ModuleBase::GlobalFunc::ZEROS(Rtau[i][j] + beg, len);
				}
			}
			
			for(int i=0; i<GlobalV::NSPIN; i++)
			{
				ModuleBase::GlobalFunc::ZEROS(tau_save2[i] + beg, len);
				for(int j=0; j<dstep; j++)
				{
					ModuleBase::GlobalFunc::ZEROS( dRtau[i][j] + beg, len );
					ModuleBase::GlobalFunc::ZEROS( dtau[i][j] + beg, len);
				}
			}			
		}
		
		ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, dstep, 512, beg, len);
		ModuleBase::GlobalFunc::ZEROS(dRR + beg, len);
		ModuleBase::GlobalFunc::ZEROS(alpha + beg, len);
	};

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
            	Rrho[is][i] = new double[this->rhopw->nrxx];
        	}
		}
    	ModuleBase::Memory::record("ChgMix::Rrho", sizeof(double) * GlobalV::NSPIN*rstep*this->rhopw->nrxx);

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
			rho_save2[is] = new double[this->rhopw->nrxx];

			for (int i=0; i<dstep; i++)
			{
				dRrho[is][i] = new double[this->rhopw->nrxx];	
				drho[is][i] = new double[this->rhopw->nrxx];
			}
		}
    	ModuleBase::Memory::record("ChgMix::dRrho", sizeof(double) * GlobalV::NSPIN*dstep*this->rhopw->nrxx);
    	ModuleBase::Memory::record("ChgMix::drho", sizeof(double) * GlobalV::NSPIN*dstep*this->rhopw->nrxx);
    	ModuleBase::Memory::record("ChgMix::rho_save2", sizeof(double) * GlobalV::NSPIN*this->rhopw->nrxx);

		if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
		{
			this->Rtau = new double**[GlobalV::NSPIN];
			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				Rtau[is] = new double*[rstep];
				for (int i=0; i<rstep; i++)
				{
					Rtau[is][i] = new double[this->rhopw->nrxx];
				}
			}
			ModuleBase::Memory::record("ChgMix::Rtau", sizeof(double) * GlobalV::NSPIN*rstep*this->rhopw->nrxx);

			this->dRtau = new double**[GlobalV::NSPIN];
			this->dtau = new double**[GlobalV::NSPIN];
			this->tau_save2 = new double*[GlobalV::NSPIN];

			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				dRtau[is] = new double*[dstep];
				dtau[is] = new double*[dstep];
				tau_save2[is] = new double[this->rhopw->nrxx];

				for (int i=0; i<dstep; i++)
				{
					dRtau[is][i] = new double[this->rhopw->nrxx];	
					dtau[is][i] = new double[this->rhopw->nrxx];
				}
			}
			ModuleBase::Memory::record("ChgMix::dRtau", sizeof(double) * GlobalV::NSPIN*dstep*this->rhopw->nrxx);
			ModuleBase::Memory::record("ChgMix::dtau", sizeof(double) * GlobalV::NSPIN*dstep*this->rhopw->nrxx);
			ModuleBase::Memory::record("ChgMix::tau_save2", sizeof(double) * GlobalV::NSPIN*this->rhopw->nrxx);			
		}

		ModuleBase::GlobalFunc::NOTE("Allocate Abar = <dRrho_j | dRrho_i >, dimension = dstep.");
		this->Abar.create(dstep, dstep);
		ModuleBase::Memory::record("ChgMix::Abar", sizeof(double) * dstep*dstep);

		// (4) allocate dRR = <delta R|R>
		ModuleBase::GlobalFunc::NOTE("Allocate dRR = < dR | R >, dimension = dstep");
		this->dRR = new double[dstep];

		// (5) allocate alpha
		ModuleBase::GlobalFunc::NOTE("Allocate alpha, dimension = dstep");
		this->alpha = new double[dstep];

		// (6) zeros all arrays
		ModuleBase::OMP_PARALLEL(zeros_kernel);

		this->initp = true;
    }

	// mohan add 2010-07-16
	if(this->new_e_iteration)
	{
		ModuleBase::OMP_PARALLEL(zeros_kernel);
	}
	return;
}

void Charge_Mixing::deallocate_Pulay()
{
    if (!this->initp) return;
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

	// dimension of rho_save2(GlobalV::NSPIN, this->rhopw->nrxx)
	for (int is=0; is<GlobalV::NSPIN; is++)
	{
		delete[] rho_save2[is];
	}
	delete[] rho_save2;	

	if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
	{
		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			for (int i=0; i<rstep; i++)
			{
				delete[] Rtau[is][i];
			}
			delete[] Rtau[is];
		}
		delete[] Rtau;

		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			for (int i=0; i<dstep; i++)
			{
				delete[] dRtau[is][i];
				delete[] dtau[is][i];
			}
			delete[] dRtau[is];
			delete[] dtau[is];
		}
		delete[] dRtau;
		delete[] dtau;

		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] tau_save2[is];
		}
		delete[] tau_save2;	
	}
    this->initp = false;
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
				A(i,j) += calculate_residual_norm(this->rhopw->nrxx, this->dRrho[is][j], this->dRrho[is][i] );
				if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
				{
					A(i,j) += calculate_residual_norm(this->rhopw->nrxx, this->dRtau[is][j], this->dRtau[is][i] );
				}
				A(j,i) = A(i,j);
			}
		}
	}
	Parallel_Reduce::reduce_double_pool(A.c, A.nr * A.nc);
	return;
}

#include "module_base/complexmatrix.h"
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
			this->dRR[i] += calculate_residual_norm(this->rhopw->nrxx, this->dRrho[is][i], this->Rrho[is][m]);
			if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
			{
				this->dRR[i] += calculate_residual_norm(this->rhopw->nrxx, this->dRtau[is][i], this->Rtau[is][m]);
			}
		}
	}
	Parallel_Reduce::reduce_double_pool(dRR, dstep);

	return;
}

// use dstep to genearte Abar(dstep, dstep)
void Charge_Mixing::generate_alpha()
{
	for(int i=0; i<dstep; i++)
	{
		this->alpha[i] = 0;
		for(int j=0; j<dstep; j++)
		{
			// Abar is a symmetry matrix.
			this->alpha[i] -= this->Abar(j,i) * this->dRR[j];
		}	
	}

	return;
}

void Charge_Mixing::generate_new_rho(const int &is, const int &m, Charge* chr)
{	
	double mixp = this->mixing_beta;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int ir=0; ir<this->rhopw->nrxx; ir++)
	{
		double rhonew = chr->rho_save[is][ir] + mixp * this->Rrho[is][m][ir];
		for(int i=0; i<dstep; i++)
		{
			rhonew += this->alpha[i] * ( this->drho[is][i][ir] + mixp * this->dRrho[is][i][ir] );
		}
		chr->rho[is][ir] = rhonew;
	}

	if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
	{
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(int ir=0; ir<this->rhopw->nrxx; ir++)
		{
			double rhonew = chr->kin_r_save[is][ir] + mixp * this->Rtau[is][m][ir];
			for(int i=0; i<dstep; i++)
			{
				rhonew += this->alpha[i] * ( this->dtau[is][i][ir] + mixp * this->dRtau[is][i][ir] );
			}
			chr->kin_r[is][ir] = rhonew;
		}
	}
	return;
}

void Charge_Mixing::generate_residual_vector(double *residual, const double* rho_out, const double* rho_in)const
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
	for (int ir=0; ir<this->rhopw->nrxx; ir++)
	{
		residual[ir]= rho_out[ir] - rho_in[ir];
	}
	return;
}


// calculate "dR^{i} = R^{i+1} - R^{i}"
// calculate "drho^{i} = rho^{i+1} - rho^{i}"
void Charge_Mixing::generate_datas(const int &irstep, const int &idstep, const int &totstep, Charge* chr)
{
	//===============================================
	// calculate the important "Rrho".
	// init the "Rrho = rho - rho_save"
	// which Rrho to be update now? answer: irstep
	//===============================================

	ModuleBase::GlobalFunc::NOTE("Generate Residual std::vector from rho and rho_save.");
    for (int is=0; is<GlobalV::NSPIN; is++)
    {
		this->generate_residual_vector( this->Rrho[is][irstep], chr->rho[is], chr->rho_save[is]);
		if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
		{
			this->generate_residual_vector( this->Rtau[is][irstep], chr->kin_r[is], chr->kin_r_save[is]);
		}

		// Note: there is no kerker modification for tau because I'm not sure
		// if we should have it. If necessary we can try it in the future.

		if(this->mixing_gg0 > 0.0)
		{
			std::complex<double> *kerpulay = new std::complex<double>[this->rhopw->npw];
			double* kerpulayR = new double[this->rhopw->nrxx];
			
			this->rhopw->real2recip(Rrho[is][irstep], kerpulay);

			const double fac = this->mixing_gg0;
			const double gg0 = std::pow(fac * 0.529177 /GlobalC::ucell.tpiba, 2);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 128)
#endif
			for(int ig=0; ig<this->rhopw->npw; ig++)
			{
				double gg = this->rhopw->gg[ig];
				double filter_g = std::max(gg / (gg + gg0), 0.1);
				kerpulay[ig] = (1 - filter_g) * kerpulay[ig];
			}

			this->rhopw->recip2real(kerpulay, kerpulayR);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
			for(int ir=0; ir<this->rhopw->nrxx; ir++)
			{
				Rrho[is][irstep][ir] = Rrho[is][irstep][ir] - kerpulayR[ir];
			}

			delete[] kerpulay;
			delete[] kerpulayR;
		}
    }

	ModuleBase::OMP_PARALLEL([&](int num_threads, int thread_id)
	{
		int irbeg, irend;
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, this->rhopw->nrxx, irbeg, irend);
		irend = irbeg + irend;
		if(totstep==0)
		{
			// don't need to calculate 'dRrho' and 'drho'.

			// mohan fix the bug 2010/03/26.
			// save 'rho_save2' in order to calculate drho in
			// the next iteration.
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				for (int ir=irbeg; ir<irend; ir++)
				{
					// save 'rho_save2' in order to calculate drho in
					this->rho_save2[is][ir] = chr->rho_save[is][ir];
				}
			}

			if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
			{
				for(int is=0; is<GlobalV::NSPIN; is++)
				{
					for (int ir=irbeg; ir<irend; ir++)
					{
						// save 'tau_save2' in order to calculate drho in
						this->tau_save2[is][ir] = chr->kin_r_save[is][ir];
					}
				}		
			}
		}
		else if(totstep>0)
		{	
			// which dRrho to be update now? answer: dstep=istep-1;
			// irstep range: [0, rstep)
			const int nowR = irstep;
			int lastR = irstep - 1;

			if (thread_id == 0)
			{
				if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "now irstep", nowR);
				if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "last irstep",lastR); 
			}

			if(lastR < 0) lastR += rstep;
			for (int is=0; is<GlobalV::NSPIN; is++)
			{
				for (int ir=irbeg; ir<irend; ir++)
				{
					this->dRrho[is][idstep][ir] = this->Rrho[is][nowR][ir] - this->Rrho[is][lastR][ir];
					this->drho[is][idstep][ir] = chr->rho_save[is][ir] - this->rho_save2[is][ir];
					// save 'rho_save2' in order to calculate drho in
					this->rho_save2[is][ir] = chr->rho_save[is][ir];
				}
			}

			if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
			{
				for (int is=0; is<GlobalV::NSPIN; is++)
				{
					for (int ir=irbeg; ir<irend; ir++)
					{
						this->dRtau[is][idstep][ir] = this->Rtau[is][nowR][ir] - this->Rtau[is][lastR][ir];
						this->dtau[is][idstep][ir] = chr->kin_r_save[is][ir] - this->tau_save2[is][ir];
						// save 'tau_save2' in order to calculate drho in
						this->tau_save2[is][ir] = chr->kin_r_save[is][ir];
					}
				}			
			}
		}
	});

	ModuleBase::GlobalFunc::NOTE("Calculate drho = rho_{in}^{i+1} - rho_{in}^{i}");
	return;
}

