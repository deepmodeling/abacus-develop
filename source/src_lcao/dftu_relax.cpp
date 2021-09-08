//==========================================================
// Author:Xin Qu
// DATE : 2019-12-10
//==========================================================
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <sstream>
#include <complex>

#include "dftu_relax.h"
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
}
namespace ModuleDFTU{
DFTU_RELAX::DFTU_RELAX(){}

DFTU_RELAX::~DFTU_RELAX(){}

void DFTU_RELAX::force_stress()
{
	ModuleBase::TITLE("DFTU_RELAX", "force_stress");

	if(GlobalV::FORCE)
	{
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			for(int dim=0; dim<3; dim++)
			{
				this->force_dftu.at(iat).at(dim) = 0.0;
			}
		}
	}

	if(GlobalV::STRESS)
	{
		for(int dim=0; dim<3; dim++)
		{
			this->stress_dftu.at(dim).at(0) = 0.0;
			this->stress_dftu.at(dim).at(1) = 0.0;
			this->stress_dftu.at(dim).at(2) = 0.0;
		}
  }

  if(GlobalV::GAMMA_ONLY_LOCAL)
  {
    const char transN = 'N', transT = 'T';
	  const int  one_int = 1;
	  const double alpha = 1.0, beta = 0.0;

    std::vector<double> rho_VU(GlobalC::ParaO.nloc);

    for(int ik=0; ik<GlobalC::kv.nks; ik++)
	  {
	  	const int spin = GlobalC::kv.isk[ik];

      double* VU = new double [GlobalC::ParaO.nloc];
      this->cal_VU_pot_mat_real(spin, false, VU);
      pdgemm_(&transT, &transN,
					    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
					    &alpha, 
					    GlobalC::LOC.wfc_dm_2d.dm_gamma[spin].c, &one_int, &one_int, GlobalC::ParaO.desc, 
					    VU, &one_int, &one_int, GlobalC::ParaO.desc,
					    &beta,
					    &rho_VU[0], &one_int, &one_int, GlobalC::ParaO.desc);

      delete [] VU;

      if(GlobalV::FORCE) this->cal_force_gamma(&rho_VU[0]);
      if(GlobalV::STRESS) this->cal_stress_gamma(&rho_VU[0]);
    }//ik
  }
  else
  {
    const char transN = 'N', transT = 'T';
	  const int  one_int = 1;
	  const std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);

    std::vector<std::complex<double>> rho_VU(GlobalC::ParaO.nloc);

    for(int ik=0; ik<GlobalC::kv.nks; ik++)
	  {
	  	const int spin = GlobalC::kv.isk[ik];

      std::complex<double>* VU = new std::complex<double> [GlobalC::ParaO.nloc];
      this->cal_VU_pot_mat_complex(spin, false, VU);
      pzgemm_(&transT, &transN,
					    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
					    &alpha, 
					    GlobalC::LOC.wfc_dm_2d.dm_k[ik].c, &one_int, &one_int, GlobalC::ParaO.desc, 
					    VU, &one_int, &one_int, GlobalC::ParaO.desc,
					    &beta,
					    &rho_VU[0], &one_int, &one_int, GlobalC::ParaO.desc);

      delete [] VU;

      if(GlobalV::FORCE) this->cal_force_k(ik, &rho_VU[0]);
      if(GlobalV::STRESS) cal_stress_k(ik, &rho_VU[0]);
    }//ik
  }

	if(GlobalV::FORCE)
	{
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			std::vector<double> tmp = this->force_dftu[iat];
			MPI_Allreduce(&tmp[0], &this->force_dftu[iat][0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		}
	}

  if(GlobalV::STRESS)
  {
    for(int dim1=0; dim1<3; dim1++)
	  {
	  	for(int dim2=dim1; dim2<3; dim2++)
	  	{
        double val_tmp = stress_dftu[dim1][dim2];
        MPI_Allreduce(&val_tmp, &stress_dftu[dim1][dim2], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
    }

    for(int i=0; i<3; i++)
	  	for(int j=0; j<3; j++)
	  		if(i>j) this->stress_dftu[i][j] = stress_dftu[j][i];

	  for(int i=0;i<3;i++)
	  	for(int j=0;j<3;j++)
	  		this->stress_dftu[i][j] *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
  }

	return;
}

void DFTU_RELAX::cal_force_k(const int ik, const std::complex<double>* rho_VU)
{
	ModuleBase::TITLE("DFTU_RELAX", "cal_force_k");

	const char transN = 'N';
	const int  one_int = 1;
  const std::complex<double> zero(0.0,0.0), one(1.0,0.0);

	std::vector<std::complex<double>> dm_VU_dSm(GlobalC::ParaO.nloc);
  std::vector<std::complex<double>> dSm_k(GlobalC::ParaO.nloc);

	for(int dim=0; dim<3; dim++)
	{
    this->fold_dSm_k(ik, dim, &dSm_k[0]);
    
		pzgemm_(&transN, &transN,
			&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
			&one, 
			rho_VU, &one_int, &one_int, GlobalC::ParaO.desc, 
			&dSm_k[0], &one_int, &one_int, GlobalC::ParaO.desc,
			&zero,
			&dm_VU_dSm[0], &one_int, &one_int, GlobalC::ParaO.desc);

    for(int ir=0; ir<GlobalC::ParaO.nrow; ir++)
		{
			const int iwt1 = GlobalC::ParaO.MatrixInfo.row_set[ir];
			const int iat1 = GlobalC::ucell.iwt2iat[iwt1];

			for(int ic=0; ic<GlobalC::ParaO.ncol; ic++)
			{
				const int iwt2 = GlobalC::ParaO.MatrixInfo.col_set[ic];
				const int irc = ic*GlobalC::ParaO.nrow + ir;

				if(iwt1==iwt2) this->force_dftu[iat1][dim] += 2.0*dm_VU_dSm[irc].real();

			}//end ic
		}//end ir

	}//end dim

	return;
}

void DFTU_RELAX::cal_stress_k(const int ik, const std::complex<double>* rho_VU)
{
	ModuleBase::TITLE("DFTU_RELAX", "cal_stress_k");

	const char transN = 'N';
	const int  one_int = 1;
	const std::complex<double> minus_half(-0.5,0.0), zero(0.0,0.0), one(1.0,0.0);

  std::vector<std::complex<double>> dm_VU_sover(GlobalC::ParaO.nloc);
  std::vector<std::complex<double>> dSR_k(GlobalC::ParaO.nloc);

	for(int dim1=0; dim1<3; dim1++)
	{
		for(int dim2=dim1; dim2<3; dim2++)
		{
			this->fold_dSR_k(ik, dim1, dim2, &dSR_k[0]);

			pzgemm_(&transN, &transN,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
				&minus_half, 
				rho_VU, &one_int, &one_int, GlobalC::ParaO.desc, 
				&dSR_k[0], &one_int, &one_int, GlobalC::ParaO.desc,
				&zero,
				&dm_VU_sover[0], &one_int, &one_int, GlobalC::ParaO.desc);

			for(int ir=0; ir<GlobalC::ParaO.nrow; ir++)
			{
				const int iwt1 = GlobalC::ParaO.MatrixInfo.row_set[ir];
				for(int ic=0; ic<GlobalC::ParaO.ncol; ic++)
				{
					const int iwt2 = GlobalC::ParaO.MatrixInfo.col_set[ic];
					const int irc = ic*GlobalC::ParaO.nrow + ir;

					if(iwt1==iwt2) stress_dftu[dim1][dim2] += 2.0*dm_VU_sover[irc].real();
				}//end ic
			}//end ir

		}//end dim2	
	}//end dim1
	
	return;
}

void DFTU_RELAX::cal_force_gamma(const double* rho_VU)
{
	ModuleBase::TITLE("DFTU_RELAX", "cal_force_gamma");

	const char transN = 'N';
	const int  one_int = 1;
	const double one = 1.0, zero = 0.0, minus_one=-1.0;
	
	std::vector<double> dm_VU_dSm(GlobalC::ParaO.nloc);
	
	for(int dim=0; dim<3; dim++)
	{
    double* tmp_ptr;
    if(dim==0) tmp_ptr = GlobalC::LM.DSloc_x;
    else if(dim==1) tmp_ptr = GlobalC::LM.DSloc_y;
    else if(dim==2) tmp_ptr = GlobalC::LM.DSloc_z;

		pdgemm_(&transN, &transN,
			&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
			&one, 
			rho_VU, &one_int, &one_int, GlobalC::ParaO.desc, 
			tmp_ptr, &one_int, &one_int, GlobalC::ParaO.desc,
			&zero,
			&dm_VU_dSm[0], &one_int, &one_int, GlobalC::ParaO.desc);

    for(int ir=0; ir<GlobalC::ParaO.nrow; ir++)
		{
			const int iwt1 = GlobalC::ParaO.MatrixInfo.row_set[ir];
			const int iat1 = GlobalC::ucell.iwt2iat[iwt1];

			for(int ic=0; ic<GlobalC::ParaO.ncol; ic++)
			{
				const int iwt2 = GlobalC::ParaO.MatrixInfo.col_set[ic];
				const int irc = ic*GlobalC::ParaO.nrow + ir;

				if(iwt1==iwt2) this->force_dftu[iat1][dim] += 2.0*dm_VU_dSm[irc];

			}//end ic
		}//end ir
	}// end dim

	return;
}

void DFTU_RELAX::cal_stress_gamma(const double* rho_VU)
{
	ModuleBase::TITLE("DFTU_RELAX", "cal_stress_gamma");

	const char transN = 'N';
	const int  one_int = 1;
	const double zero = 0.0, minus_half=-0.5, one=1.0;

  std::vector<double> dSR_gamma(GlobalC::ParaO.nloc);
  std::vector<double> dm_VU_sover(GlobalC::ParaO.nloc);

	for(int dim1=0; dim1<3; dim1++)
	{
		for(int dim2=dim1; dim2<3; dim2++)
		{	
      this->fold_dSR_gamma(dim1, dim2, &dSR_gamma[0]);

			pdgemm_(&transN, &transN,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
				&minus_half, 
				rho_VU, &one_int, &one_int, GlobalC::ParaO.desc, 
				&dSR_gamma[0], &one_int, &one_int, GlobalC::ParaO.desc,
				&zero,
				&dm_VU_sover[0], &one_int, &one_int, GlobalC::ParaO.desc);

			for(int ir=0; ir<GlobalC::ParaO.nrow; ir++)
			{
				const int iwt1 = GlobalC::ParaO.MatrixInfo.row_set[ir];

				for(int ic=0; ic<GlobalC::ParaO.ncol; ic++)
				{
					const int iwt2 = GlobalC::ParaO.MatrixInfo.col_set[ic];
					const int irc = ic*GlobalC::ParaO.nrow + ir;

					if(iwt1==iwt2) stress_dftu[dim1][dim2] += 2.0*dm_VU_sover[irc];
				}//end ic
			}//end ir

		}//end dim2
	}//end dim1

	return;
}

double DFTU_RELAX::get_onebody_eff_pot
(
	const int T, const int iat,
	const int L, const int N, const int spin, 
	const int m0, const int m1,
  const int type, const bool newlocale 
)
{
	ModuleBase::TITLE("DFTU_RELAX","get_onebody_eff_pot");

	double VU = 0.0;

	//if(!Yukawa && N!=0) return 0.0;

	switch(type)
	{
	case 1:  //rotationally invarient formalism and FLL double counting
		/*
		const int spin_oppsite = 1 - spin;
		int nelec_tot = 0;
		int nelec_spin = 0;

		for(int is=0; is<2; is++)
		{
			for(int i=0; i<2*L+1; i++)
			{
				if(newlocale)
				{
					nelec_tot +=  locale.at(iat).at(L).at(N).at(is)(i,i);
					if(is==spin) nelec_spin +=  locale.at(iat).at(L).at(N).at(is)(i,i);
				}
				else
				{
					nelec_tot +=  locale_save.at(iat).at(L).at(N).at(is)(i,i);
					if(is==spin) nelec_spin +=  locale_save.at(iat).at(L).at(N).at(is)(i,i);
				}
				
			}
		}

		for(int m2=0; m2<2*L+1; m2++)
		{
			for(int m3=0; m3<2*L+1; m3++)
			{
				int M0_prime = m2*(2*L+1) + m3;

				int M1 = m0*(2*L+1) + m3;
				int M1_prime = m2*(2*L+1) + m1;
				
				if(newlocale)
				{
					VU += Vsc.at(T).at(N)(M0,M0_prime)*locale.at(iat).at(L).at(N).at(spin_oppsite)(m2,m3) 
						+ (Vsc.at(T).at(N)(M0,M0_prime)- Vsc.at(T).at(N)(M1,M1_prime))*locale.at(iat).at(L).at(N).at(spin)(m2,m3);
				}
				else
				{
					VU += Vsc.at(T).at(N)(M0, M0_prime)*locale_save.at(iat).at(L).at(N).at(spin_oppsite)(m2,m3) 
						+ (Vsc.at(T).at(N)(M0, M0_prime)- Vsc.at(T).at(N)(M1,M1_prime))*locale_save.at(iat).at(L).at(N).at(spin)(m2,m3);
				}									
			}
		}

		//if(Yukawa) VU += -U_Yukawa.at(T).at(N)*(nelec_tot-0.5) + J_Yukawa.at(T).at(N)*(nelec_spin-0.5);
		//else VU += -U[T]*(nelec_tot-0.5) + J[T]*(nelec_spin-0.5);
		VU += -U[T]*(nelec_tot-0.5) + J[T]*(nelec_spin-0.5);
		*/

		break;

	case 2:  //rotationally invarient formalism and AMF double counting
		
		break;

	case 3:	//simplified formalism and FLL double counting
		if(newlocale)
		{
			if(Yukawa)
			{				
				if(m0==m1) VU = (this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*( 0.5 - this->locale.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*this->locale.at(iat).at(L).at(N).at(spin)(m0,m1);
			}
			else
			{
				if(m0==m1) VU = (this->U[T]-this->J[T])*( 0.5 - this->locale.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U[T]-this->J[T])*this->locale.at(iat).at(L).at(N).at(spin)(m0,m1);
			}	
		}
		else
		{	
			if(Yukawa)
			{				
				if(m0==m1) VU = (this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*( 0.5 - this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U_Yukawa.at(T).at(L).at(N)-this->J_Yukawa.at(T).at(L).at(N))*this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1);
			}
			else
			{
				if(m0==m1) VU = (this->U[T]-this->J[T])*( 0.5 - this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1) );		
				else VU = -(this->U[T]-this->J[T])*this->locale_save.at(iat).at(L).at(N).at(spin)(m0,m1);
			}
		}
		
		break;

	case 4: //simplified formalism and AMF double counting
		
		break;		
	}

	return VU;
}

void DFTU_RELAX::cal_VU_pot_mat_complex(const int spin, const bool newlocale, std::complex<double>* VU)
{
  ModuleBase::TITLE("DFTU_RELAX","cal_VU_pot_mat_complex"); 
	// ModuleBase::timer::tick("DFTU","folding_overlap_matrix");
  ModuleBase::GlobalFunc::ZEROS(VU, GlobalC::ParaO.nloc);

  for(int it=0; it<GlobalC::ucell.ntype; ++it)
	{
    if(INPUT.orbital_corr[it]==-1) continue;
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
      const int iat = GlobalC::ucell.itia2iat(it, ia);
			for(int L=0; L<=GlobalC::ucell.atoms[it].nwl; L++)
			{			
        if(L!=INPUT.orbital_corr[it] ) continue;

				for(int n=0; n<GlobalC::ucell.atoms[it].l_nchi[L]; n++)
				{
					// if(Yukawa)
			    // {
			    	// if(L1<INPUT.orbital_corr[T1] || L2<INPUT.orbital_corr[T2]) continue;
			    // }
			    // else
			    // {
			    	if(n!=0) continue;
			    // }
          for(int m1=0; m1<2*L+1; m1++)
          {
            for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
		  		  {
		  			  const int mu = GlobalC::ParaO.trace_loc_row[this->iatlnmipol2iwt[iat][L][n][m1][ipol1]];
              if(mu<0) continue;

              for(int m2=0; m2<2*L+1; m2++)
              {
                for(int ipol2=0; ipol2<GlobalV::NPOL; ipol2++)
                {
                  const int nu = GlobalC::ParaO.trace_loc_col[this->iatlnmipol2iwt[iat][L][n][m2][ipol2]];
                  if(nu<0) continue;

                  int m1_all = m1 + (2*L+1)*ipol1;
			            int m2_all = m2 + (2*L+1)*ipol2;
                  
                  double val = get_onebody_eff_pot(it, iat, L, n, spin, m1_all, m2_all, cal_type, newlocale);
			            VU[nu*GlobalC::ParaO.nrow + mu] = std::complex<double>(val, 0.0);
                }//ipol2
              }//m2
            }//ipol1
          }//m1
				}//n
			}//l
		}//ia	
	}//it

  return;
}

void DFTU_RELAX::cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU)
{
  ModuleBase::TITLE("DFTU_RELAX","cal_VU_pot_mat_real"); 
	// ModuleBase::timer::tick("DFTU","folding_overlap_matrix");
  ModuleBase::GlobalFunc::ZEROS(VU, GlobalC::ParaO.nloc);

  for(int it=0; it<GlobalC::ucell.ntype; ++it)
	{
    if(INPUT.orbital_corr[it]==-1) continue;
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
      const int iat = GlobalC::ucell.itia2iat(it, ia);
			for(int L=0; L<=GlobalC::ucell.atoms[it].nwl; L++)
			{			
        if(L!=INPUT.orbital_corr[it] ) continue;

				for(int n=0; n<GlobalC::ucell.atoms[it].l_nchi[L]; n++)
				{
					// if(Yukawa)
			    // {
			    	// if(L1<INPUT.orbital_corr[T1] || L2<INPUT.orbital_corr[T2]) continue;
			    // }
			    // else
			    // {
			    	if(n!=0) continue;
			    // }
          for(int m1=0; m1<2*L+1; m1++)
          {
            for(int ipol1=0; ipol1<GlobalV::NPOL; ipol1++)
		  		  {
		  			  const int mu = GlobalC::ParaO.trace_loc_row[this->iatlnmipol2iwt[iat][L][n][m1][ipol1]];
              if(mu<0) continue;
              for(int m2=0; m2<2*L+1; m2++)
              {
                for(int ipol2=0; ipol2<GlobalV::NPOL; ipol2++)
                {
                  const int nu = GlobalC::ParaO.trace_loc_col[this->iatlnmipol2iwt[iat][L][n][m2][ipol2]];
                  if(nu<0) continue;

                  int m1_all = m1 + (2*L+1)*ipol1;
			            int m2_all = m2 + (2*L+1)*ipol2;
                  
                  VU[nu*GlobalC::ParaO.nrow + mu] = this->get_onebody_eff_pot(it, iat, L, n, spin, m1_all, m2_all, cal_type, newlocale);

                }//ipol2
              }//m2
            }//ipol1
          }//m1
				}//n
			}//l
		}//ia	
	}//it

  return;
}

void DFTU_RELAX::fold_dSR_gamma(const int dim1, const int dim2, double* dSR_gamma)
{
  ModuleBase::TITLE("DFTU_RELAX","fold_dSR_gamma");

  ModuleBase::GlobalFunc::ZEROS(dSR_gamma, GlobalC::ParaO.nloc);

  double* dS_ptr;
  if(dim1==0) dS_ptr =  GlobalC::LM.DSloc_x;
  else if(dim1==1) dS_ptr =  GlobalC::LM.DSloc_y;
  else if(dim1==2) dS_ptr =  GlobalC::LM.DSloc_z;

  int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
  for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
  {
	  Atom* atom1 = &GlobalC::ucell.atoms[T1];
    for(int I1=0; I1<atom1->na; ++I1)
    {
		  tau1 = atom1->tau[I1];
      const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);    
      GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
      for(int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1; ++ad)
      {
        const int T2 = GlobalC::GridD.getType(ad);
			  const int I2 = GlobalC::GridD.getNatom(ad);
        const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
			  Atom* atom2 = &GlobalC::ucell.atoms[T2];
			  tau2 = GlobalC::GridD.getAdjacentTau(ad);
			  dtau = tau2 - tau1;
			  double distance = dtau.norm() * GlobalC::ucell.lat0;
			  double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
        bool adj = false;
			  if(distance < rcut) adj = true;
			  else if(distance >= rcut)
			  {
			  	for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
			  	{
			  		const int T0 = GlobalC::GridD.getType(ad0); 
			  		const int I0 = GlobalC::GridD.getNatom(ad0); 
			  		const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
			  		const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
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

			  if(adj)
			  {
			  	for(int jj=0; jj<atom1->nw*GlobalV::NPOL; ++jj)
			  	{
            const int jj0 = jj/GlobalV::NPOL;
            const int iw1_all = start1 + jj0; 
            const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
				    if(mu<0)continue;

			  		for(int kk=0; kk<atom2->nw*GlobalV::NPOL; ++kk)
			  		{
              const int kk0 = kk/GlobalV::NPOL;
              const int iw2_all = start2 + kk0;
					    const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
					    if(nu<0)continue;

              dSR_gamma[nu*GlobalC::ParaO.nrow + mu] += dS_ptr[nnr]*GlobalC::LM.DH_r[nnr*3+dim2];

			  			++nnr;
			  		}// kk
			    }// jj
			  }// adj
		  }// ad
	  }// I1
	}// T1

	return;
}

void DFTU_RELAX::fold_dSm_k(const int ik, const int dim, std::complex<double>* dSm_k)
{
  ModuleBase::TITLE("DFTU_RELAX","fold_dSm_k");

  ModuleBase::GlobalFunc::ZEROS(dSm_k, GlobalC::ParaO.nloc);

  double* dSm_ptr;
  if(dim==0) dSm_ptr = GlobalC::LM.DSloc_Rx;
  else if(dim==1) dSm_ptr = GlobalC::LM.DSloc_Ry;
  else if(dim==2) dSm_ptr = GlobalC::LM.DSloc_Rz;

  int nnr = 0;
  ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
  for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
  {
	  Atom* atom1 = &GlobalC::ucell.atoms[T1];
    for(int I1=0; I1<atom1->na; ++I1)
    {
		  tau1 = atom1->tau[I1];
      const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);    

      GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
      for(int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1; ++ad)
      {
        	const int T2 = GlobalC::GridD.getType(ad);
			  const int I2 = GlobalC::GridD.getNatom(ad);
        const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

			  Atom* atom2 = &GlobalC::ucell.atoms[T2];

			  tau2 = GlobalC::GridD.getAdjacentTau(ad);
			  dtau = tau2 - tau1;

			  double distance = dtau.norm() * GlobalC::ucell.lat0;
			  double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

        bool adj = false;
			  if(distance < rcut) adj = true;
			  else if(distance >= rcut)
			  {
			  	for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
			  	{
			  		const int T0 = GlobalC::GridD.getType(ad0); 
			  		const int I0 = GlobalC::GridD.getNatom(ad0); 
			  		const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
			  		const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

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

			  if(adj)
			  {
			  	for(int jj=0; jj<atom1->nw*GlobalV::NPOL; ++jj)
			  	{
					const int jj0 = jj/GlobalV::NPOL;
					const int iw1_all = start1 + jj0; 
					const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
				    if(mu<0) continue;

			  		for(int kk=0; kk<atom2->nw*GlobalV::NPOL; ++kk)
			  		{
              			const int kk0 = kk/GlobalV::NPOL;
              			const int iw2_all = start2 + kk0;
					    const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
					    if(nu<0) continue;

			  			ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z); 
			  			const double arg = ( GlobalC::kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
			  			const std::complex<double> kphase( cos(arg),  sin(arg) );

              			dSm_k[nu*GlobalC::ParaO.nrow + mu] += dSm_ptr[nnr]*kphase;

			  			++nnr;
			  		}// kk
			    }// jj
			  }// adj

		  }// ad
	  }// I1
	}// T1

	return;
}

void DFTU_RELAX::fold_dSR_k(const int ik, const int dim1, const int dim2, std::complex<double>* dSR_k)
{
  ModuleBase::TITLE("DFTU_RELAX","fold_dSR_k");

  ModuleBase::GlobalFunc::ZEROS(dSR_k, GlobalC::ParaO.nloc);

  double* dSm_ptr;
  if(dim1==0) dSm_ptr = GlobalC::LM.DSloc_Rx;
  else if(dim1==1) dSm_ptr = GlobalC::LM.DSloc_Ry;
  else if(dim1==2) dSm_ptr = GlobalC::LM.DSloc_Rz;

  int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
  for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
  {
	  Atom* atom1 = &GlobalC::ucell.atoms[T1];
    for(int I1=0; I1<atom1->na; ++I1)
    {
		  tau1 = atom1->tau[I1];
      const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);    

      GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
      for(int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1; ++ad)
      {
        const int T2 = GlobalC::GridD.getType(ad);
			  const int I2 = GlobalC::GridD.getNatom(ad);
        const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

			  Atom* atom2 = &GlobalC::ucell.atoms[T2];

			  tau2 = GlobalC::GridD.getAdjacentTau(ad);
			  dtau = tau2 - tau1;

			  double distance = dtau.norm() * GlobalC::ucell.lat0;
			  double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

        bool adj = false;
			  if(distance < rcut) adj = true;
			  else if(distance >= rcut)
			  {
			  	for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
			  	{
			  		const int T0 = GlobalC::GridD.getType(ad0); 
			  		const int I0 = GlobalC::GridD.getNatom(ad0); 
			  		const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
			  		const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

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

			  if(adj)
			  {
			  	for(int jj=0; jj<atom1->nw*GlobalV::NPOL; ++jj)
			  	{
            const int jj0 = jj/GlobalV::NPOL;
            const int iw1_all = start1 + jj0; 
            const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
				    if(mu<0)continue;

			  		for(int kk=0; kk<atom2->nw*GlobalV::NPOL; ++kk)
			  		{
              const int kk0 = kk/GlobalV::NPOL;
              const int iw2_all = start2 + kk0;
					    const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
					    if(nu<0)continue;
              	
			  			ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z); 
			  			const double arg = ( GlobalC::kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
			  			const std::complex<double> kphase( cos(arg),  sin(arg) );

			  			dSR_k[nu*GlobalC::ParaO.nrow + mu] += dSm_ptr[nnr]*GlobalC::LM.DH_r[nnr*3+dim2]*kphase;														

			  			++nnr;
			  		}// kk
			    }// jj
			  }// adj
        
		  }// ad
	  }// I1
	}// T1

	return;
}

}