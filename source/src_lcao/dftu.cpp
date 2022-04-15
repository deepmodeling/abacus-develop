//==========================================================
// Author:Xin Qu
// DATE : 2019-12-10
//==========================================================
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <sstream>
#include <complex>
#include "../module_base/memory.h"

#include "dftu.h"
#include "../module_base/constants.h"
#include "../src_pw/global.h"
#include "global_fp.h"
#include "../module_base/global_function.h"
#include "../module_base/inverse_matrix.h"
#include "LOOP_ions.h"
#include "LCAO_matrix.h"
#include "../src_pw/magnetism.h"
#include "../module_orbital/ORB_gen_tables.h"
#include "../src_pw/charge.h"
#include "../module_base/timer.h"

extern "C"
{
  void pzgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const std::complex<double> *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const std::complex<double> *beta,
		std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
  
  void pdgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC);

  void pztranc_(
    const int *M, const int *N,
    const std::complex<double> *alpha,
    const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
    const std::complex<double> *beta,
    std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);

  void pdtran_(
    const int *M, const int *N,
    const double *alpha,
    const double *A, const int *IA, const int *JA, const int *DESCA,
    const double *beta,
    double *C, const int *IC, const int *JC, const int *DESCC);
}

namespace GlobalC
{
ModuleDFTU::DFTU dftu;
}

namespace ModuleDFTU{
DFTU::DFTU(){}

DFTU::~DFTU(){}

void DFTU::init(
	UnitCell_pseudo &cell, // unitcell class
    LCAO_Matrix &lm
    )
{
    ModuleBase::TITLE("DFTU", "init");

    this->LM = &lm;
    
    // needs reconstructions in future
	// global parameters, need to be removed in future
	const int npol = GlobalV::NPOL; // number of polarization directions
	const int nlocal = GlobalV::NLOCAL; // number of total local orbitals
	const int nks = GlobalC::kv.nks; // number of k-points
	const int nspin = GlobalV::NSPIN; // number of spins
	
	

	if(dftu_type==1 && double_counting==1) cal_type = 1;
	else if(dftu_type==1 && double_counting==2) cal_type = 2;
	else if(dftu_type==2 && double_counting==1) cal_type = 3;
	else if(dftu_type==2 && double_counting==2) cal_type = 4;
	else ModuleBase::WARNING_QUIT("DFT+U", "Wrong parameter");
		
	this->EU = 0.0;

	if(GlobalV::CAL_FORCE)
	{
		this->force_dftu.resize(cell.nat);
		for(int ia=0; ia<cell.nat; ia++)
			this->force_dftu.at(ia).resize(3, 0.0);
	}

	if(GlobalV::CAL_STRESS)
	{
		this->stress_dftu.resize(3);
		for(int dim=0; dim<3; dim++)
			this->stress_dftu.at(dim).resize(3, 0.0);
	}

	this->locale.resize(cell.nat);
	this->locale_save.resize(cell.nat);

	this->iatlnmipol2iwt.resize(cell.nat);
	this->iat2it.resize(cell.nat);
	this->iwt2it.resize(nlocal);
	this->iwt2l.resize(nlocal);
	this->iwt2n.resize(nlocal);
	this->iwt2m.resize(nlocal);
	this->iwt2ipol.resize(nlocal);	

	//initialize 
	for(int i=0; i<nlocal; i++)
	{
		this->iwt2it.at(i) = -1;
		this->iwt2l.at(i) = -1;
		this->iwt2n.at(i) = -1;
		this->iwt2m.at(i) = -1;
		this->iwt2ipol.at(i) = -1;
	}

	int num_locale=0;
	//it:index of type of atom
	for(int it=0; it<cell.ntype; ++it)
	{				
		for(int ia=0; ia<cell.atoms[it].na; ia++)
		{
			//ia:index of atoms of this type
			//determine the size of locale
			const int iat = cell.itia2iat(it, ia);
			this->iat2it.at(iat) = it;

			locale.at(iat).resize(cell.atoms[it].nwl+1);
			locale_save.at(iat).resize(cell.atoms[it].nwl+1);

			for(int l=0; l<=cell.atoms[it].nwl; l++)
			{			
				const int N = cell.atoms[it].l_nchi[l];

				locale.at(iat).at(l).resize(N);
				locale_save.at(iat).at(l).resize(N);

				for(int n=0; n<N; n++)
				{
					if(nspin==1 || nspin==2)
					{
						locale.at(iat).at(l).at(n).resize(2);
						locale_save.at(iat).at(l).at(n).resize(2);

						locale.at(iat).at(l).at(n).at(0).create(2*l+1, 2*l+1);
						locale.at(iat).at(l).at(n).at(1).create(2*l+1, 2*l+1);

						locale_save.at(iat).at(l).at(n).at(0).create(2*l+1, 2*l+1);
						locale_save.at(iat).at(l).at(n).at(1).create(2*l+1, 2*l+1);
						num_locale+=(2*l+1)*(2*l+1)*2;
					}
					else if(nspin==4) //SOC
					{
						locale.at(iat).at(l).at(n).resize(1);
						locale_save.at(iat).at(l).at(n).resize(1);

						locale.at(iat).at(l).at(n).at(0).create((2*l+1)*npol, (2*l+1)*npol);
						locale_save.at(iat).at(l).at(n).at(0).create((2*l+1)*npol, (2*l+1)*npol);
						num_locale+=(2*l+1)*(2*l+1)*npol*npol;
					}												
				}
			}

			//initialize the arrry iatlnm2iwt[iat][l][n][m]
			this->iatlnmipol2iwt.at(iat).resize(cell.atoms[it].nwl+1);
			for(int L=0; L<=cell.atoms[it].nwl; L++)
			{
				this->iatlnmipol2iwt.at(iat).at(L).resize(cell.atoms[it].l_nchi[L]);

				for(int n=0; n<cell.atoms[it].l_nchi[L]; n++)
				{
					this->iatlnmipol2iwt.at(iat).at(L).at(n).resize(2*L+1);
					
					for(int m=0; m<2*L+1; m++)
					{
						this->iatlnmipol2iwt.at(iat).at(L).at(n).at(m).resize(npol);
					}
				}
			}

			for(int iw=0; iw<cell.atoms[it].nw*npol; iw++)
			{
				int iw0 = iw/npol;
				int ipol = iw%npol;
				int iwt = cell.itiaiw2iwt(it, ia, iw);
				int l = cell.atoms[it].iw2l[iw0];
				int n = cell.atoms[it].iw2n[iw0];
				int m = cell.atoms[it].iw2m[iw0];
								
				this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m).at(ipol) = iwt;
				this->iwt2it.at(iwt) = it;
				this->iwt2l.at(iwt) = l;
				this->iwt2n.at(iwt) = n;
				this->iwt2m.at(iwt) = m;
				this->iwt2ipol.at(iwt) = ipol;
			}
		}	
	}

	
 	if(Yukawa)
 	{
		this->Fk.resize(cell.ntype);
		
		this->U_Yukawa.resize(cell.ntype);
		this->J_Yukawa.resize(cell.ntype);	

		for(int it=0; it<cell.ntype; it++)
		{
			const int NL = cell.atoms[it].nwl + 1;

			this->Fk.at(it).resize(NL);		
			this->U_Yukawa.at(it).resize(NL);
			this->J_Yukawa.at(it).resize(NL);	

			for(int l=0; l<NL; l++)
			{
				int N = cell.atoms[it].l_nchi[l];

				this->Fk.at(it).at(l).resize(N);
				for(int n=0; n<N; n++)
				{
					this->Fk.at(it).at(l).at(n).resize(l+1, 0.0);
				}	

				this->U_Yukawa.at(it).at(l).resize(N, 0.0);
				this->J_Yukawa.at(it).at(l).resize(N, 0.0);

				// if(l>=orbital_corr[it] && orbital_corr[it]!=-1)
				// {
					// this->cal_slater_Fk(l, it);
					// this->cal_slater_UJ(it, l);
				// }
			}			 	
		}
 	}

	if(GlobalV::CALCULATION=="nscf" || GlobalC::pot.init_chg == "file")
	{
		std::stringstream sst; 
		sst << GlobalV::global_out_dir << "onsite.dm"; 
		this->read_occup_m( sst.str() );
		this->local_occup_bcast();		
	}
	else
	{
		if(omc) 
		{
			std::stringstream sst; 
			sst << "initial_onsite.dm"; 
			this->read_occup_m( sst.str() );
			this->local_occup_bcast();
		}
	}
	
	ModuleBase::Memory::record("DFTU","locale",num_locale,"double");

  //GlobalV::ofs_running << "GlobalC::dftu.cpp "<< __LINE__ << std::endl;
    return;
}

void DFTU::cal_occup_m_k(const int iter,  std::vector<ModuleBase::ComplexMatrix> &dm_k)
{
	ModuleBase::TITLE("DFTU", "cal_occup_m_k");
	ModuleBase::timer::tick("DFTU", "cal_occup_m_k");

	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{
		if(orbital_corr[T]==-1) continue;

		//const int LC = orbital_corr[T];
		
		for(int I=0; I<GlobalC::ucell.atoms[T].na; I++)
		{
		 	const int iat = GlobalC::ucell.itia2iat(T, I);
			
			for(int l=0; l<GlobalC::ucell.atoms[T].nwl+1; l++)
			{
				const int N = GlobalC::ucell.atoms[T].l_nchi[l];

				for(int n=0; n<N; n++)
				{						
					if(GlobalV::NSPIN==4)
					{
						locale_save[iat][l][n][0] = locale[iat][l][n][0];
            			locale[iat][l][n][0].zero_out();
					}
					else if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
					{
						locale_save[iat][l][n][0] = locale[iat][l][n][0];
						locale_save[iat][l][n][1] = locale[iat][l][n][1];
						locale[iat][l][n][0].zero_out();
						locale[iat][l][n][1].zero_out();
					}
				}
			}			
		}
	}

	//=================Part 1======================
	//call PBLAS routine to calculate the product of the S and density matrix
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);

	std::vector<std::complex<double>> srho(this->LM->ParaV->nloc);
    std::vector<std::complex<double>> Sk(this->LM->ParaV->nloc);
	
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		// srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
        this->folding_overlap_matrix(ik, &Sk[0], *this->LM);

	#ifdef __MPI
		pzgemm_(&transN, &transT,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
				&alpha, 
				&Sk[0], &one_int, &one_int, this->LM->ParaV->desc, 
				dm_k.at(ik).c, &one_int, &one_int, this->LM->ParaV->desc,
				&beta, 
				&srho[0], &one_int, &one_int, this->LM->ParaV->desc);
	#endif

    	const int spin = GlobalC::kv.isk[ik];
    	for(int it=0; it<GlobalC::ucell.ntype; it++)
	  	{
			const int NL = GlobalC::ucell.atoms[it].nwl + 1;
			const int LC = orbital_corr[it];
	
			if(LC == -1) continue;

		  	for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		  	{		
		  	const int iat = GlobalC::ucell.itia2iat(it, ia);

		  	for(int l=0; l<NL; l++)
		  	{				
		  		// if(Yukawa)
		  		// {
		  			// if(l<orbital_corr[it]) continue;
		  		// }
		  		// else
		  		// {
		  			// if(l!=orbital_corr[it]) continue;
		  		// }
		  		if(l!=orbital_corr[it]) continue;

		  		const int N = GlobalC::ucell.atoms[it].l_nchi[l];
  
		  		for(int n=0; n<N; n++)
		  		{
		  		 	// if(!Yukawa && n!=0) continue;
		  			if(n!=0) continue;

		  			//Calculate the local occupation number matrix			
		  			for(int m0=0; m0<2*l+1; m0++)
		  			{
		  				for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
		  				{
		  					const int iwt0 = this->iatlnmipol2iwt[iat][l][n][m0][ipol0];
		  					const int mu = this->LM->ParaV->trace_loc_row[iwt0];
		  					const int mu_prime = this->LM->ParaV->trace_loc_col[iwt0];

		  					for(int m1=0; m1<2*l+1; m1++)
		  					{
		  						for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
		  						{									
		  							const int iwt1 = this->iatlnmipol2iwt[iat][l][n][m1][ipol1];
		  							const int nu = this->LM->ParaV->trace_loc_col[iwt1];
		  							const int nu_prime = this->LM->ParaV->trace_loc_row[iwt1];

		  							const int irc = nu*this->LM->ParaV->nrow + mu;
		  							const int irc_prime = mu_prime*this->LM->ParaV->nrow + nu_prime;

		  							const int m0_all = m0 + ipol0*(2*l+1);
		  							const int m1_all = m1 + ipol1*(2*l+1);

		  							if( (nu>=0) && (mu>=0) )
		  								locale[iat][l][n][spin](m0_all, m1_all) += (srho[irc]).real()/4.0;									

		  							if( (nu_prime>=0) && (mu_prime>=0) )
		  								locale[iat][l][n][spin](m0_all, m1_all) += (std::conj(srho[irc_prime])).real()/4.0;
		  						}//ipol1										
		  					}//m1
		  				}//ipol0
		  			}//m0
		  			}//end n
		  		// this->print(it, iat, l, N, iter);
		  		}//end l
		  	}//end ia
	  	}//end it
	}//ik

  	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
	  	const int NL = GlobalC::ucell.atoms[it].nwl + 1;
	  	const int LC = orbital_corr[it];
  
	  	if(LC == -1) continue;

		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{		
			const int iat = GlobalC::ucell.itia2iat(it, ia);

			for(int l=0; l<NL; l++)
			{				
				// if(Yukawa)
				// {
					// if(l<orbital_corr[it]) continue;
				// }
				// else
				// {
					// if(l!=orbital_corr[it]) continue;
				// }
				if(l!=orbital_corr[it]) continue;

	  			const int N = GlobalC::ucell.atoms[it].l_nchi[l];

				for(int n=0; n<N; n++)
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;
					// set the local occupation mumber matrix of spin up and down zeros

				#ifdef __MPI
					if(GlobalV::NSPIN==1 || GlobalV::NSPIN==4)
					{
						ModuleBase::matrix temp(locale[iat][l][n][0]);
						MPI_Allreduce( &temp(0,0), &locale[iat][l][n][0](0,0), (2*l+1)*GlobalV::NPOL*(2*l+1)*GlobalV::NPOL,
												MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
					else if(GlobalV::NSPIN==2)
					{
						ModuleBase::matrix temp0(locale[iat][l][n][0]);
						MPI_Allreduce( &temp0(0,0), &locale[iat][l][n][0](0,0), (2*l+1)*(2*l+1),
												MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

						ModuleBase::matrix temp1(locale[iat][l][n][1]);
						MPI_Allreduce( &temp1(0,0), &locale[iat][l][n][1](0,0), (2*l+1)*(2*l+1),
												MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
				#endif

					// for the case spin independent calculation
					switch(GlobalV::NSPIN)
					{
						case 1:
						locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
						locale[iat][l][n][0] *= 0.5;
						locale[iat][l][n][1] += locale[iat][l][n][0];
						break;

						case 2:
						for(int is=0; is<GlobalV::NSPIN; is++)
							locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
						break;

						case 4: //SOC
						locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
						break;

						default:
						std::cout << "Not supported NSPIN parameter" << std::endl;
						exit(0);			
					}
				}//end n
	  		// this->print(it, iat, l, N, iter);
	  		}//end l
	  	}//end ia
	}//end it

	//GlobalV::ofs_running << "GlobalC::dftu.cpp "<< __LINE__  << std::endl;

	ModuleBase::timer::tick("DFTU", "cal_occup_m_k");
	return;
}

void DFTU::cal_occup_m_gamma(const int iter,  std::vector<ModuleBase::matrix> &dm_gamma)
{
	ModuleBase::TITLE("DFTU", "cal_occup_m_gamma");	
	ModuleBase::timer::tick("DFTU", "cal_occup_m_gamma");
	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{
		if(orbital_corr[T]==-1) continue;
		const int L = orbital_corr[T];
		
		for(int I=0; I<GlobalC::ucell.atoms[T].na; I++)
		{
			const int iat = GlobalC::ucell.itia2iat(T, I);
			
			for(int l=L; l<GlobalC::ucell.atoms[T].nwl+1; l++)
			{
				const int N = GlobalC::ucell.atoms[T].l_nchi[l];

				for(int n=0; n<N; n++)
				{									
					locale_save[iat][l][n][0] = locale[iat][l][n][0];
					locale_save[iat][l][n][1] = locale[iat][l][n][1];	
					locale[iat][l][n][0].zero_out();
					locale[iat][l][n][1].zero_out();				
				}
			}			
		}
	}
	
	//=================Part 1======================
	//call PBLAS routine to calculate the product of the S and density matrix
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;

	std::vector<double> srho(this->LM->ParaV->nloc);
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		// srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_gamma(iw,nu)
	
	#ifdef __MPI
		pdgemm_(&transN, &transT,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
				&alpha, 
				this->LM->Sloc.data(), &one_int, &one_int, this->LM->ParaV->desc, 
				dm_gamma.at(is).c, &one_int, &one_int, this->LM->ParaV->desc,
				&beta,
				&srho[0], &one_int, &one_int, this->LM->ParaV->desc);
	#endif

    for(int it=0; it<GlobalC::ucell.ntype; it++)
	  {
	  	const int NL = GlobalC::ucell.atoms[it].nwl + 1;
	  	if(orbital_corr[it] == -1) continue;
	  	for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
	  	{
	  		const int iat = GlobalC::ucell.itia2iat(it, ia);

	  		for(int l=0; l<NL; l++)
	  		{				
	  			// if(Yukawa)
	  			// {
	  				// if(l<orbital_corr[it]) continue;
	  			// }
	  			// else
	  			// {
	  				// if(l!=orbital_corr[it]) continue;
	  			// }
	  			if(l!=orbital_corr[it]) continue;

	  			const int N = GlobalC::ucell.atoms[it].l_nchi[l];

	  			for(int n=0; n<N; n++)
	  			{
	  			 	// if(!Yukawa && n!=0) continue;
	  				if(n!=0) continue;

	  				//Calculate the local occupation number matrix			
	  				for(int m0=0; m0<2*l+1; m0++)
	  				{	
	  					for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
	  					{
	  						const int iwt0 = this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m0).at(ipol0);
	  						const int mu = this->LM->ParaV->trace_loc_row[iwt0];
	  						const int mu_prime = this->LM->ParaV->trace_loc_col[iwt0];

	  						for(int m1=0; m1<2*l+1; m1++)
	  						{	
	  							for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
	  							{											
	  								const int iwt1 = this->iatlnmipol2iwt.at(iat).at(l).at(n).at(m1).at(ipol1);
	  								const int nu = this->LM->ParaV->trace_loc_col[iwt1];
	  								const int nu_prime = this->LM->ParaV->trace_loc_row[iwt1];

	  								const int irc = nu*this->LM->ParaV->nrow + mu;
	  								const int irc_prime = mu_prime*this->LM->ParaV->nrow + nu_prime;

	  								if( (nu>=0) && (mu>=0) )
	  								{																																																
	  									int m0_all = m0 + (2*l+1)*ipol0;
	  									int m1_all = m0 + (2*l+1)*ipol1;

	  									locale[iat][l][n][is](m0,m1) += srho[irc]/4.0;														
	  								}

	  								if( (nu_prime>=0) && (mu_prime>=0) )
	  								{
	  									int m0_all = m0 + (2*l+1)*ipol0;
	  									int m1_all = m0 + (2*l+1)*ipol1;
  
	  									locale[iat][l][n][is](m0,m1) += srho[irc_prime]/4.0;
	  								}
	  							}			
	  						}
	  					}
	  				}

	  				ModuleBase::matrix temp(locale[iat][l][n][is]);
					
				#ifdef __MPI
	  				MPI_Allreduce( &temp(0,0), &locale[iat][l][n][is](0,0), (2*l+1)*GlobalV::NPOL*(2*l+1)*GlobalV::NPOL,
	  								      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				#endif

	  				// for the case spin independent calculation
	  				switch(GlobalV::NSPIN)
	  				{
	  				  case 1:	
	  				  	locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                locale[iat][l][n][0] *= 0.5;
                locale[iat][l][n][1] += locale[iat][l][n][0];
	  				  	break;

	  				  case 2:
	  				  	locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
	  				  	break;

	  				  default:
	  				  	std::cout << "Not supported NSPIN parameter" << std::endl;
	  				  	exit(0);			
	  				}

	  			}//end for(n)
	  			//this->print(it, iat, 2*L+1, N);
	  		}//L
	  	}//ia
	  }//it

	}//is

	//test the sum rule of srho
	/*
	double elec_tot = 0.0;
	double Nele = 0.0; 
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int ir=0; ir<this->LM->ParaV->nrow; ir++)
		{
			for(int ic=0; ic<this->LM->ParaV->ncol; ic++)
			{
				int row = this->LM->ParaV->MatrixInfo.row_set[ir];
				int col = this->LM->ParaV->MatrixInfo.col_set[ic];

				if(row==col) elec_tot += srho.at(is).at(ic*this->LM->ParaV->nrow + ir);
			}
		}
	}	
	MPI_Allreduce(&elec_tot, &Nele, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if(GlobalV::MY_RANK==0)
	{
		std::ofstream ofs_elec("nelec.dat", ios_base::app);
		ofs_elec << "Total number of electrons of the system " << Nele << std::endl;
	}
	*/

	//ofs_running << "GlobalC::dftu.cpp "<< __LINE__  << std::endl;
	ModuleBase::timer::tick("DFTU", "cal_occup_m_gamma");
	return;
}

void DFTU::write_occup_m(const std::string &fn)
{
	ModuleBase::TITLE("DFTU", "write_occup_m");

	if(GlobalV::MY_RANK!=0) return;

	std::ofstream ofdftu;
	ofdftu.open(fn.c_str());
    if (!ofdftu)
    {
    	std::cout << "DFTU::write_occup_m. Can't create file onsite.dm!" << std::endl;
		exit(0);
    }

	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{
		if(orbital_corr[T]==-1) continue;
		const int NL = GlobalC::ucell.atoms[T].nwl+1;
		const int LC = orbital_corr[T];

		for(int I=0; I<GlobalC::ucell.atoms[T].na; I++)
		{
			const int iat = GlobalC::ucell.itia2iat(T, I);
			ofdftu << "atoms" << "  " << iat << std::endl;

			for(int l=0; l<NL; l++)
			{
				// if(Yukawa)
				// {
					// if(l<orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=orbital_corr[T]) continue;
				// }
				if(l!=orbital_corr[T]) continue;
				
				const int N = GlobalC::ucell.atoms[T].l_nchi[l];
				ofdftu << "L" << "  " << l << std::endl;

				for(int n=0; n<N; n++)		
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

					ofdftu << "zeta" << "  " << n << std::endl;
					
					if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
					{
						for(int is=0; is<2; is++)
						{
							ofdftu << "spin" << "  " << is << std::endl;
							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int m1=0; m1<2*l+1; m1++)
								{
									ofdftu << fixed << std::setw(12) << std::setprecision(8) << locale.at(iat).at(l).at(n).at(is)(m0, m1);
								}
								ofdftu << std::endl;
							}
						}
					}
					else if(GlobalV::NSPIN==4) //SOC
					{
						for(int m0=0; m0<2*l+1; m0++)
						{
							for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
							{
								const int m0_all = m0 + (2*l+1)*ipol0;

								for(int m1=0; m1<2*l+1; m1++)
								{
									for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
									{
										int m1_all = m1 + (2*l+1)*ipol1;
										ofdftu << fixed << std::setw(12) << std::setprecision(8) << locale.at(iat).at(l).at(n).at(0)(m0_all, m1_all);
									}
								}
								ofdftu << std::endl;
							}							
						}
					}					
				}//n
			}//l
		}//I
	}//T

	return;
}

void DFTU::read_occup_m(const std::string &fn)
{
	ModuleBase::TITLE("DFTU", "read_occup_m");

	if(GlobalV::MY_RANK!=0) return;

	std::ifstream ifdftu(fn.c_str(), ios::in);	

	if (!ifdftu) 
	{
		if(GlobalV::CALCULATION=="nscf") 
		{
			std::cout << "DFTU::read_occup_m. Can not find the file onsite.dm . Please do scf calculation first" << std::endl;
		}
		else
		{
			if(omc) 
			{
				std::cout << "DFTU::read_occup_m. Can not find the file initial_onsite.dm . Please check your initial_onsite.dm" << std::endl;
			}
		}
		exit(0);
	}

	ifdftu.clear();
    ifdftu.seekg(0);

	char word[10];

	int T, iat, spin, L, zeta;

	ifdftu.rdstate();

    while(ifdftu.good())
    {
        ifdftu >> word;
		if(ifdftu.eof()) break;
		
		if(strcmp("atoms", word) == 0)
        {
            ifdftu >> iat;
			ifdftu.ignore(150, '\n');

			T= this->iat2it.at(iat);
			const int NL = GlobalC::ucell.atoms[T].nwl + 1;
			const int LC = orbital_corr[T];
	
			for(int l=0; l<NL; l++)			
			{
				// if(Yukawa)
				// {
					// if(l<orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=orbital_corr[T]) continue;
				// }
				if(l!=orbital_corr[T]) continue;

				ifdftu >> word;

				if (strcmp("L", word) == 0)
				{
					ifdftu >> L;
					ifdftu.ignore(150, '\n');

					const int N = GlobalC::ucell.atoms[T].l_nchi[L];
					for(int n=0; n<N; n++)
					{
						// if(!Yukawa && n!=0) continue;
						if(n!=0) continue;

						ifdftu >> word;
						if(strcmp("zeta", word) == 0)						
						{														
							ifdftu >> zeta;
							ifdftu.ignore(150, '\n');
							
							if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
							{
								for(int is=0; is<2; is++)
								{
									ifdftu >> word;
									if(strcmp("spin", word) == 0)
									{								
										ifdftu >> spin;
										ifdftu.ignore(150, '\n');

										double value = 0.0;
										for(int m0=0; m0<2*L+1; m0++)
										{
											for(int m1=0; m1<2*L+1; m1++)
											{							
												ifdftu >> value;
												locale.at(iat).at(L).at(zeta).at(spin)(m0, m1) = value;						 	
											}
											ifdftu.ignore(150, '\n');											
										}
									}
									else
									{								
										std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << std::endl;
										exit(0);
									}	
								}
							}
							else if(GlobalV::NSPIN==4) //SOC
							{
								double value = 0.0;
								for(int m0=0; m0<2*L+1; m0++)
								{
									for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
									{
										const int m0_all = m0 + (2*L+1)*ipol0;

										for(int m1=0; m1<2*L+1; m1++)
										{		
											for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)	
											{		
												int m1_all = m1 + (2*L+1)*ipol1;	
												ifdftu >> value;
												locale.at(iat).at(L).at(zeta).at(0)(m0_all, m1_all) = value;		
											}				 	
										}
										ifdftu.ignore(150, '\n');
									}											
								}
							}																																		
						}
						else
						{						
							std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << std::endl;
							exit(0);
						}
					}
				}
				else
				{				
					std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << std::endl;
					exit(0);
				}								
			}			
        }
		else
		{		
			std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM DFTU FILE" << std::endl;
			exit(0);
		}
		
		ifdftu.rdstate();

		if (ifdftu.eof() != 0)
  		{
			break;
  		}
	}
   
	return;
}

void DFTU::local_occup_bcast()
{
	ModuleBase::TITLE("DFTU", "local_occup_bcast");

	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{
		if(orbital_corr[T]==-1) continue;

		for(int I=0; I<GlobalC::ucell.atoms[T].na; I++)
		{
			const int iat = GlobalC::ucell.itia2iat(T,I);
			const int L = orbital_corr[T];

			for(int l=0; l<=GlobalC::ucell.atoms[T].nwl; l++)
			{
				// if(Yukawa)
				// {
					// if(l<orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=orbital_corr[T]) continue;
				// }
				if(l!=orbital_corr[T]) continue;

				for(int n=0; n<GlobalC::ucell.atoms[T].l_nchi[l]; n++)
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

				#ifdef __MPI				
					if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
					{
						for(int spin=0; spin<2; spin++)
						{
						   for(int m0=0; m0<2*l+1; m0++)
						   {
							   for(int m1=0; m1<2*l+1; m1++)
							   {
									MPI_Bcast(&locale.at(iat).at(l).at(n).at(spin)(m0, m1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
							   }
						   }
						}
					}
					else if(GlobalV::NSPIN==4) //SOC
					{
						for(int m0=0; m0<2*L+1; m0++)
						{
							for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
							{
								const int m0_all = m0 + (2*L+1)*ipol0;

								for(int m1=0; m1<2*L+1; m1++)
								{		
									for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)	
									{		
										int m1_all = m1 + (2*L+1)*ipol1;	
										
										MPI_Bcast(&locale.at(iat).at(l).at(n).at(0)(m0_all, m1_all), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);		
									}				 	
								}
							}											
						}
					}		
				#endif			
				}
			}
		}
	}	    
	return;
}

void DFTU::cal_energy_correction(const int istep)
{
	ModuleBase::TITLE("DFTU", "cal_energy_correction");
	ModuleBase::timer::tick("DFTU", "cal_energy_correction");
	if((GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax") && (!omc) && istep==0 && this->iter_dftu==1 && GlobalC::pot.init_chg != "file")
	{
		ModuleBase::timer::tick("DFTU", "cal_energy_correction"); 
		return;
	}
	this->EU = 0.0;
	double EU_dc = 0.0;
	
	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{
		const int NL = GlobalC::ucell.atoms[T].nwl + 1;
		const int LC = orbital_corr[T]; 
		for(int I=0; I<GlobalC::ucell.atoms[T].na; I++)
		{
			if(LC == -1) continue;

			const int iat = GlobalC::ucell.itia2iat(T, I);
			const int L = orbital_corr[T];

			for(int l=0; l<NL; l++)
			{
				// if(Yukawa)
				// {
					// if(l<orbital_corr[T]) continue;
				// }	
				// else if(l!=orbital_corr[T]) continue;
				if(l!=orbital_corr[T]) continue;

				const int N = GlobalC::ucell.atoms[T].l_nchi[l];

				const int m_tot = 2*l+1;

				//part 1: calculate the DFT+U energy correction
				for(int n=0; n<N; n++)
				{
				 	// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

					if(cal_type==1) //rotationally invarient formalism and FLL double counting; Not available at present
					{

					}
					else if(cal_type==2) //rotationally invarient formalism and AMF double counting; Not available at present
					{

					}
					else if(cal_type==3) //simplified formalism and FLL double counting
					{	
						if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
						{
							for(int spin=0; spin<2; spin++)
							{
								double nm_trace = 0.0;
								double nm2_trace = 0.0;

								for(int m0=0; m0<2*l+1; m0++)
								{
								 	nm_trace += this->locale.at(iat).at(l).at(n).at(spin)(m0, m0);
									for(int m1=0; m1<2*l+1; m1++)
									{
										nm2_trace += this->locale.at(iat).at(l).at(n).at(spin)(m0,m1)*this->locale.at(iat).at(l).at(n).at(spin)(m1,m0);
									}
								}
								if(Yukawa) this->EU += 0.5*(this->U_Yukawa.at(T).at(l).at(n) - this->J_Yukawa.at(T).at(l).at(n))*(nm_trace - nm2_trace);
								else this->EU += 0.5*(this->U[T] - this->J[T])*(nm_trace - nm2_trace);	
							}
						}
						else if(GlobalV::NSPIN==4) //SOC
						{
							double nm_trace = 0.0;
							double nm2_trace = 0.0;

							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
								{
							 		const int m0_all = m0 + (2*l+1)*ipol0;
									nm_trace += this->locale.at(iat).at(l).at(n).at(0)(m0_all, m0_all);

									for(int m1=0; m1<2*l+1; m1++)
									{
										for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
										{
											int m1_all = m1 + (2*l+1)*ipol1;

											nm2_trace += this->locale.at(iat).at(l).at(n).at(0)(m0_all,m1_all)*this->locale.at(iat).at(l).at(n).at(0)(m1_all, m0_all);
										}
									}
								}
							}
							if(Yukawa) this->EU += 0.5*(this->U_Yukawa.at(T).at(l).at(n) - this->J_Yukawa.at(T).at(l).at(n))*(nm_trace - nm2_trace);
							else this->EU += 0.5*(this->U[T] - this->J[T])*(nm_trace - nm2_trace);
						}
						
					}
					else if(cal_type==4) //simplified formalism and AMF double counting; ot available at present
					{

					}

					//calculate the double counting term included in eband													
					for(int m1=0; m1<2*l+1; m1++)
					{
						for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
						{
							const int m1_all = m1 + ipol1*(2*l+1);
							for(int m2=0; m2<2*l+1; m2++)
							{	
								for(int ipol2=0; ipol2<GlobalV::NPOL; ipol2++)
								{
									const int m2_all = m2 + ipol2*(2*l+1);

									if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
									{
										for(int is=0; is<2; is++)
										{
											double VU = 0.0;
											VU = get_onebody_eff_pot(T, iat, l, n, is, m1_all, m2_all, cal_type, 0);
											EU_dc += VU*this->locale.at(iat).at(l).at(n).at(is)(m1_all,m2_all);
										}
									}	
									else if(GlobalV::NSPIN==4) //SOC
									{										
										double VU = 0.0;
										VU = get_onebody_eff_pot(T, iat, l, n, 0, m1_all, m2_all, cal_type, 0);
										EU_dc += VU*this->locale.at(iat).at(l).at(n).at(0)(m1_all,m2_all);										
									}
								}					
							}	
						}
					}						
				}//end n
			}//end L
		}//end I
	}//end T
	
	/*
	if(GlobalV::MY_RANK==0)
	{
		std::ofstream of_eu("energy_correction", ios_base::app);
		double e=EU*ModuleBase::Ry_to_eV;
		of_eu << "ITERATION STEP " << std::setw(3) << this->iter_dftu <<"      " << std::setw(15) << e << "ev" << std::endl;
	}
	*/
	
	//substract the double counting EU_dc included in band energy eband
	this->EU -= EU_dc;

	/*
	if(GlobalV::MY_RANK==0)
	{
		std::ofstream of_en("energy_correction.dat",ios_base::app);
		double val = this->EU*ModuleBase::Ry_to_eV;

		of_en << fixed << std::setprecision(8) << val << "eV" << std::endl;
	}
	*/
	ModuleBase::timer::tick("DFTU", "cal_energy_correction");
	return;
}

void DFTU::cal_eff_pot_mat_complex(const int ik, const int istep, std::complex<double>* eff_pot)
{
	ModuleBase::TITLE("DFTU", "cal_eff_pot_mat");
	ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
 	if((GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax") && (!omc) && istep==0 && this->iter_dftu==1 && GlobalC::pot.init_chg != "file")
	{
		ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
		return;
	} 
	
	int spin = GlobalC::kv.isk[ik];

	ModuleBase::GlobalFunc::ZEROS(eff_pot, this->LM->ParaV->nloc);

	//GlobalV::ofs_running << "GlobalC::dftu.cpp "<< __LINE__  << std::endl;
	//=============================================================
	//   PART2: call pblas to calculate effective potential matrix
	//=============================================================
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
  	const std::complex<double> zero(0.0,0.0), half(0.5,0.0), one(1.0,0.0);

	std::vector<std::complex<double>> VU(this->LM->ParaV->nloc);
  	this->cal_VU_pot_mat_complex(spin, true, &VU[0]);

#ifdef __MPI
	pzgemm_(&transN, &transN,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&half, 
		ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->LM->ParaV->desc,
		this->LM->Sloc2.data(), &one_int, &one_int, this->LM->ParaV->desc,
		&zero,
		eff_pot, &one_int, &one_int, this->LM->ParaV->desc);
#endif

  	for(int irc=0; irc<this->LM->ParaV->nloc; irc++)
    VU[irc] = eff_pot[irc];
  
  	//pztranc(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

#ifdef __MPI
  	pztranc_(&GlobalV::NLOCAL, &GlobalV::NLOCAL, 
		&one, 
		&VU[0], &one_int, &one_int, this->LM->ParaV->desc, 
		&one, 
		eff_pot, &one_int, &one_int, this->LM->ParaV->desc);
#endif

	//code for testing whther the effective potential is Hermitian
	/*
	bool pot_Hermitian = true;
	if(GAMMA_ONLY_LOCAL)
	{
		for(int i=0; i<NLOCAL; i++)
		{
			for(int j=0; j<NLOCAL; j++)
			{
				int iic = i*NLOCAL + j;
				int jjc = j*NLOCAL + i;
				double tmp = pot_eff_gamma.at(spin).at(iic) - pot_eff_gamma.at(spin).at(jjc);
				if(tmp>1.0e-9) pot_Hermitian = false;
			}
		}
	}
	else
	{
		for(int i=0; i<NLOCAL; i++)
		{
			for(int j=i; j<NLOCAL; j++)
			{
				int iic = i*NLOCAL + j;
				int jjc = j*NLOCAL + i;
				std::complex<double> tmp = pot_eff_k.at(ik).at(iic) - conj(pot_eff_k.at(ik).at(jjc));
				double tmp_norm = sqrt(std::norm(tmp));
				if(tmp_norm>1.0e-9) pot_Hermitian = false;
			}
		}
	}


	if(MY_RANK==0)
	{
		std::ofstream of_potH("Hermitian_pot.dat",ios_base::app);
		of_potH << "Hermitian  " << pot_Hermitian << std::endl;
	}
	*/

	//ofs_running << "GlobalC::dftu.cpp "<< __LINE__  << std::endl;
	ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
	return;	
}

void DFTU::cal_eff_pot_mat_real(const int ik, const int istep, double* eff_pot)
{
	ModuleBase::TITLE("DFTU", "cal_eff_pot_mat");
	ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
 	if((GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax") && (!omc) && istep==0 && this->iter_dftu==1 && GlobalC::pot.init_chg != "file")
	{
		ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
		return;
	} 
	
	int spin = GlobalC::kv.isk[ik];

	ModuleBase::GlobalFunc::ZEROS(eff_pot, this->LM->ParaV->nloc);

	//ofs_running << "GlobalC::dftu.cpp "<< __LINE__  << std::endl;
	//=============================================================
	//   PART2: call pblas to calculate effective potential matrix
	//=============================================================
	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0, half=0.5, one=1.0;

	std::vector<double> VU(this->LM->ParaV->nloc);
  	this->cal_VU_pot_mat_real(spin, 1, &VU[0]);

#ifdef __MPI
	pdgemm_(&transN, &transN,
		&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
		&half, 
		ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->LM->ParaV->desc, 
		this->LM->Sloc.data(), &one_int, &one_int, this->LM->ParaV->desc,
		&beta,
		eff_pot, &one_int, &one_int, this->LM->ParaV->desc);
#endif

  	for(int irc=0; irc<this->LM->ParaV->nloc; irc++)
    	VU[irc] = eff_pot[irc];
  
	// pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

#ifdef __MPI
	pdtran_(&GlobalV::NLOCAL, &GlobalV::NLOCAL, 
		&one, 
		&VU[0], &one_int, &one_int, this->LM->ParaV->desc, 
		&one, 
		eff_pot, &one_int, &one_int, this->LM->ParaV->desc);
#endif

	//code for testing whther the effective potential is Hermitian
	/*
	bool pot_Hermitian = true;
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		for(int i=0; i<GlobalV::NLOCAL; i++)
		{
			for(int j=0; j<GlobalV::NLOCAL; j++)
			{
				int iic = i*GlobalV::NLOCAL + j;
				int jjc = j*GlobalV::NLOCAL + i;
				double tmp = pot_eff_gamma.at(spin).at(iic) - pot_eff_gamma.at(spin).at(jjc);
				if(tmp>1.0e-9) pot_Hermitian = false;
			}
		}
	}
	else
	{
		for(int i=0; i<GlobalV::NLOCAL; i++)
		{
			for(int j=i; j<GlobalV::NLOCAL; j++)
			{
				int iic = i*GlobalV::NLOCAL + j;
				int jjc = j*GlobalV::NLOCAL + i;
				std::complex<double> tmp = pot_eff_k.at(ik).at(iic) - conj(pot_eff_k.at(ik).at(jjc));
				double tmp_norm = sqrt(std::norm(tmp));
				if(tmp_norm>1.0e-9) pot_Hermitian = false;
			}
		}
	}


	if(GlobalV::MY_RANK==0)
	{
		std::ofstream of_potH("Hermitian_pot.dat",ios_base::app);
		of_potH << "Hermitian  " << pot_Hermitian << std::endl;
	}
	*/
	
	//ofs_running << "GlobalC::dftu.cpp "<< __LINE__  << std::endl;
	ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
	return;	
}

void DFTU::output()
{
	ModuleBase::TITLE("DFTU", "output");

    if(GlobalC::CHR.out_chg)
	{
		std::stringstream sst;
		sst << GlobalV::global_out_dir << "onsite.dm";
		this->write_occup_m( sst.str() );
	}
	
	GlobalV::ofs_running << "//=========================L(S)DA+U===========================//" << std::endl;

	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{			
		const int NL = GlobalC::ucell.atoms[T].nwl + 1;

		for(int L=0; L<NL; L++)
		{
			const int N = GlobalC::ucell.atoms[T].l_nchi[L];

			if(L>=orbital_corr[T] && orbital_corr[T]!=-1)
			{
 				if(L!=orbital_corr[T]) continue;
				
				if(!Yukawa)
				{
					GlobalV::ofs_running << "atom_type=" << T << "  L=" << L << "  chi=" << 0 << "    U=" << this->U[T]*ModuleBase::Ry_to_eV 
					<< "ev    J=" << this->J[T]*ModuleBase::Ry_to_eV << "ev" << std::endl;
				}
				else
				{
					for(int n=0; n<N; n++)
					{
 						if(n!=0) continue;						
						double Ueff = (this->U_Yukawa.at(T).at(L).at(n) - this->J_Yukawa.at(T).at(L).at(n))*ModuleBase::Ry_to_eV;
						GlobalV::ofs_running << "atom_type=" << T << "  L=" << L << "  chi=" << n << "    U=" << this->U_Yukawa.at(T).at(L).at(n)*ModuleBase::Ry_to_eV << "ev    " << "J=" << this->J_Yukawa.at(T).at(L).at(n)*ModuleBase::Ry_to_eV
						<< "ev" << std::endl;
					}
				}
			}
		}
	}

	GlobalV::ofs_running << "Local occupation matrices" << std::endl;
	for(int T=0; T<GlobalC::ucell.ntype; T++)
	{
		if(orbital_corr[T]==-1) continue;
		const int NL = GlobalC::ucell.atoms[T].nwl+1;
		const int LC = orbital_corr[T];

		for(int I=0; I<GlobalC::ucell.atoms[T].na; I++)
		{
			const int iat = GlobalC::ucell.itia2iat(T, I);
			GlobalV::ofs_running << "atoms" << " " << iat ;

			for(int l=0; l<NL; l++)
			{
				// if(Yukawa)
				// {
					// if(l<orbital_corr[T]) continue;
				// }
				// else
				// {
					// if(l!=orbital_corr[T]) continue;
				// }
				if(l!=orbital_corr[T]) continue;
				
				const int N = GlobalC::ucell.atoms[T].l_nchi[l];
				GlobalV::ofs_running << "   L" << " " << l ;

				for(int n=0; n<N; n++)		
				{
					// if(!Yukawa && n!=0) continue;
					if(n!=0) continue;

					GlobalV::ofs_running << "   zeta" << " " << n << std::endl;
					
					if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
					{
						for(int is=0; is<2; is++)
						{
							GlobalV::ofs_running << "spin" << "  " << is << std::endl;
							for(int m0=0; m0<2*l+1; m0++)
							{
								for(int m1=0; m1<2*l+1; m1++)
								{
									GlobalV::ofs_running << fixed << std::setw(12) << std::setprecision(8) << locale.at(iat).at(l).at(n).at(is)(m0, m1);
								}
								GlobalV::ofs_running << std::endl;
							}
						}	
					}
					else if(GlobalV::NSPIN==4) //SOC
					{
						for(int m0=0; m0<2*l+1; m0++)
						{
							for(int ipol0=0; ipol0<GlobalV::NPOL; ipol0++)
							{
								const int m0_all = m0 + (2*l+1)*ipol0;

								for(int m1=0; m1<2*l+1; m1++)
								{
									for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
									{
										int m1_all = m1 + (2*l+1)*ipol1;
										GlobalV::ofs_running << fixed << std::setw(12) << std::setprecision(8) << locale.at(iat).at(l).at(n).at(0)(m0_all, m1_all);
									}
								}
								GlobalV::ofs_running << std::endl;
							}							
						}
					}					
				}//n
			}//l
		}//I
	}//T
	GlobalV::ofs_running << "//=======================================================//" << std::endl;
	
	return;
}

void DFTU::cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR)
{
  const char transN = 'N', transT = 'T';
  const int  one_int = 1;
  const double alpha = 1.0, beta = 0.0, one=1.0, half=0.5;

	std::vector<double> VU(this->LM->ParaV->nloc);
	this->cal_VU_pot_mat_real(ispin, 1, &VU[0]);

#ifdef __MPI
  pdgemm_(&transN, &transN,
    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
    &half, 
    ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->LM->ParaV->desc, 
    SR, &one_int, &one_int, this->LM->ParaV->desc,
    &beta,
    HR, &one_int, &one_int, this->LM->ParaV->desc);

  pdgemm_(&transN, &transN,
    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
    &half, 
    SR, &one_int, &one_int, this->LM->ParaV->desc, 
    ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->LM->ParaV->desc,
    &one,
    HR, &one_int, &one_int, this->LM->ParaV->desc);
#endif

  return;
}

void DFTU::cal_eff_pot_mat_R_complex_double(
  const int ispin, std::complex<double>* SR, std::complex<double>* HR)
{
  const char transN = 'N', transT = 'T';
  const int  one_int = 1;
  const std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);
  const std::complex<double> zero(0.0,0.0), half(0.5,0.0), one(1.0,0.0);

  std::vector<std::complex<double>> VU(this->LM->ParaV->nloc);
  this->cal_VU_pot_mat_complex(ispin, 1, &VU[0]);

#ifdef __MPI
  pzgemm_(&transN, &transN,
    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
    &half, 
    ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->LM->ParaV->desc,
    SR, &one_int, &one_int, this->LM->ParaV->desc,
    &beta,
    HR, &one_int, &one_int, this->LM->ParaV->desc);

  pzgemm_(&transN, &transN,
    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
    &half, 
    SR, &one_int, &one_int, this->LM->ParaV->desc, 
    ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->LM->ParaV->desc,
    &one,
    HR, &one_int, &one_int, this->LM->ParaV->desc);
#endif

    return;
}

void DFTU::folding_overlap_matrix(const int ik, std::complex<double>* Sk,
    LCAO_Matrix &lm)
{
  	ModuleBase::TITLE("DFTU","folding_overlap_matrix"); 
	ModuleBase::timer::tick("DFTU","folding_overlap_matrix");

  	ModuleBase::GlobalFunc::ZEROS(Sk, lm.ParaV->nloc);

	int iat = 0;
	int index = 0;
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;

	ModuleBase::Vector3<double> dtau1;
	ModuleBase::Vector3<double> dtau2;
	ModuleBase::Vector3<double> tau0;

	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GlobalC::GridD.Find_atom(tau1);
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			Atom* atom1 = &GlobalC::ucell.atoms[T1];
			const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				Atom* atom2 = &GlobalC::ucell.atoms[T2];

				tau2 = GlobalC::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				bool adj = false;

				if(distance < rcut) 
				{
					adj = true;
				}
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0); 
						const int I0 = GlobalC::GridD.getNatom(ad0); 
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							adj = true;
							break;
						}
					}
				}

				if(adj) // mohan fix bug 2011-06-26, should not be '<='
				{
					// (3) calculate the nu of atom (T2, I2)
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
					//------------------------------------------------
					// exp(k dot dR)
					// dR is the index of box in Crystal coordinates
					//------------------------------------------------
					ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z); 
					const double arg = ( GlobalC::kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
					//const double arg = ( kv.kvec_d[ik] * GlobalC::GridD.getBox(ad) ) * ModuleBase::TWO_PI;
					const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );

					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					//--------------------------------------------------
					for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
					{
						// the index of orbitals in this processor
						const int iw1_all = start + ii;
						const int mu = lm.ParaV->trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
						{
							int iw2_all = start2 + jj;
							const int nu = lm.ParaV->trace_loc_col[iw2_all];

							if(nu<0)continue;
							//const int iic = mu*lm.ParaV->ncol+nu;
							int iic;
							if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
							{
								iic=mu+nu*lm.ParaV->nrow;
							}
							else
							{
								iic=mu*lm.ParaV->ncol+nu;
							}

							//########################### EXPLAIN ###############################
							// 1. overlap matrix with k point
							// lm.SlocR = < phi_0i | phi_Rj >, where 0, R are the cell index
							// while i,j are the orbital index.

							// 2. H_fixed=T+Vnl matrix element with k point (if Vna is not used).
							// H_fixed=T+Vnl+Vna matrix element with k point (if Vna is used).
							// lm.Hloc_fixed = < phi_0i | H_fixed | phi_Rj>

							// 3. H(k) |psi(k)> = S(k) | psi(k)> 
							// Sloc2 is used to diagonalize for a give k point.
							// Hloc_fixed2 is used to diagonalize (eliminate index R).
							//###################################################################
							
							if(GlobalV::NSPIN!=4)
							{
								Sk[iic] += lm.SlocR[index] * kphase;
							}
							else
							{
								Sk[iic] += lm.SlocR_soc[index] * kphase;
							}
							++index;

						}//end jj
					}//end ii
				}
			}// end ad
			++iat;
		}// end I1
	} // end T1

	assert(index==lm.ParaV->nnr);

  	ModuleBase::timer::tick("DFTU","folding_overlap_matrix");
	return;
}

}