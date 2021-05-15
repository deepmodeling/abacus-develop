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
#include "../src_global/constants.h"
#include "../src_pw/global.h"
#include "global_fp.h"
#include "../src_global/global_function.h"
//#include "../src_global/scalapack_connector.h"
#include "../src_global/lapack_connector.h"
#include "../src_global/inverse_matrix.h"
#include "LOOP_ions.h"
#include "LCAO_matrix.h"
#include "../src_pw/magnetism.h"
#include "../module_ORB/ORB_gen_tables.h"
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

DFTU_RELAX::DFTU_RELAX(){}


DFTU_RELAX::~DFTU_RELAX(){}


void DFTU_RELAX::force_stress()
{
	TITLE("DFTU_RELAX", "force_stress");

	if(FORCE)
	{
		for(int iat=0; iat<ucell.nat; iat++)
		{
			for(int dim=0; dim<3; dim++)
			{
				this->force_dftu.at(iat).at(dim) = 0.0;
			}
		}
	}

	if(STRESS)
	{
		for(int dim=0; dim<3; dim++)
		{
			this->stress_dftu.at(dim).at(0) = 0.0;
			this->stress_dftu.at(dim).at(1) = 0.0;
			this->stress_dftu.at(dim).at(2) = 0.0;
		}
	}

    this->allocate_force_stress();

	//=======================================================
	// intialize single particel effective potential matrix
	//=======================================================
	vector<vector<complex<double>>> VU_k;
	vector<vector<double>> VU_gamma;	
	if(FORCE || STRESS)
	{
		this->folding_dSm_soverlap();

		if(GAMMA_ONLY_LOCAL)
		{
			VU_gamma.resize(NSPIN);
			for(int is=0; is<NSPIN; is++)
			{
				VU_gamma.at(is).resize(ParaO.nloc, 0.0);
			}
		}
		else
		{
			if(NSPIN==1)
			{
				VU_k.resize(1);
				VU_k.at(0).resize(ParaO.nloc, complex<double>(0.0,0.0));
			}
			else if(NSPIN==2 || NSPIN==4)
			{
				VU_k.resize(2);
				for(int is=0; is<2; is++)
				{
					VU_k.at(is).resize(ParaO.nloc, complex<double>(0.0,0.0));
				}
			}
		}

		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int T1 = this->iwt2it.at(iwt1);
			const int iat1 = this->iwt2iat.at(iwt1);
			const int L1 = this->iwt2l.at(iwt1);
			const int n1 = this->iwt2n.at(iwt1);
			const int m1 = this->iwt2m.at(iwt1);
			const int ipol1 = this->iwt2ipol.at(iwt1);			

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int T2 = this->iwt2it.at(iwt2);
				const int iat2 = this->iwt2iat.at(iwt2);
				const int L2 = this->iwt2l.at(iwt2);
				const int n2 = this->iwt2n.at(iwt2);
				const int m2 = this->iwt2m.at(iwt2);
				const int ipol2 = this->iwt2ipol.at(iwt2);				

				int irc = ic*ParaO.nrow + ir;

				if(INPUT.orbital_corr[T1]==-1 || INPUT.orbital_corr[T2]==-1) continue;

				if(iat1!=iat2) continue;

				// if(Yukawa)
				// {
					// if(L1<INPUT.orbital_corr[T1] || L2<INPUT.orbital_corr[T2]) continue;
				// }
				// else
				// {
					if(L1!=INPUT.orbital_corr[T1] || L2!=INPUT.orbital_corr[T2] || n1!=0 || n2!=0) continue;
				// }

				if(L1!=L2 || n1!=n2) continue;

				// const int m1_all = m1 + (2*L1+1)*ipol1;
				// const int m2_all = m2 + (2*L2+1)*ipol2;

				if(GAMMA_ONLY_LOCAL)
				{
					for(int is=0; is<NSPIN; is++) //no soc
					{
						double val = get_onebody_eff_pot(T1, iat1, L1, n1, is, m1, m2, cal_type, false);
					
						VU_gamma.at(is).at(irc) = val;
					}
				}
				else
				{
					if(NSPIN==1)
					{
						double val = get_onebody_eff_pot(T1, iat1, L1, n1, 0, m1, m2, cal_type, false);
						VU_k.at(0).at(irc) = complex<double>(val, 0.0);
					}
					else if(NSPIN==2)
					{
						for(int is=0; is<NSPIN; is++)
						{
							double val = get_onebody_eff_pot(T1, iat1, L1, n1, is, m1, m2, cal_type, false);
							VU_k.at(is).at(irc) = complex<double>(val, 0.0);
						}
					}
					else if(NSPIN==4)//SOC
					{
						if(ipol1==ipol2)
						{
							int is = ipol1;
							double val = get_onebody_eff_pot(T1, iat1, L1, n1, is, m1, m2, cal_type, false);
							VU_k.at(is).at(irc) = complex<double>(val, 0.0);
						}
					}
				}					
									
			}//end ic
		}//end ir
	}

	//=======================================================
	//       calculate force
	//=======================================================
	if(FORCE) 
	{
		if(GAMMA_ONLY_LOCAL) cal_force_gamma(VU_gamma);
		else  cal_force_k(VU_k);
	}

	if(STRESS)
	{
		if(GAMMA_ONLY_LOCAL) cal_stress_gamma(VU_gamma);
		else  cal_stress_k(VU_k);
	}			

	return;
}


void DFTU_RELAX::cal_force_k(const vector<vector<complex<double>>>& VU)
{
	TITLE("DFTU_RELAX", "cal_force_k");

	const char transN = 'N', transT = 'T', transC='C';
	const int  one_int = 1;
	// const double alpha = 1.0, beta = 0.0;
  const complex<double> alpha(1.0,0.0), beta(0.0,0.0);
  const complex<double> zero(0.0,0.0);
	
	vector<vector<complex<double>>> ftmp(ucell.nat);
	for(int ia=0; ia<ucell.nat; ia++)
	{
		ftmp.at(ia).resize(3, complex<double>(0.0,0.0)); //three dimension
	}

	vector<vector<complex<double>>> dm_VU_dSm(3);
	for(int dim=0; dim<3; dim++)
	{
		dm_VU_dSm.at(dim).resize(ParaO.nloc, zero);
	}
	
	for(int ik=0; ik<kv.nks; ik++)	
	{
		const int spin = kv.isk[ik];

		for(int dim=0; dim<3; dim++)
		{
			vector<complex<double>> mat_tmp(ParaO.nloc);
			vector<complex<double>> force_tmp(ParaO.nloc);

			if(dim==0) //dim=1,2 are same as dim=0
			{
				ZEROS(VECTOR_TO_PTR(mat_tmp), ParaO.nloc);

				pzgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_k.at(ik).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);
			}			

			pzgemm_(&transN, &transN,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
				this->dSm_k[ik][dim], &one_int, &one_int, ParaO.desc,
				&beta,
				VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) += force_tmp.at(irc);
			}

			//=========================================
			//   the second part
			//=========================================
			ZEROS(VECTOR_TO_PTR(force_tmp), ParaO.nloc);

			pzgemm_(&transN, &transC,
				&NLOCAL, &NLOCAL, &NLOCAL,
				&alpha, 
				this->dSm_k[ik][dim], &one_int, &one_int, ParaO.desc, 
				VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
				&beta,
				VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) -= force_tmp.at(irc);
			}
		}//end dim
	}//end ik

	for(int dim=0; dim<3; dim++)
	{
		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int iat1 = this->iwt2iat.at(iwt1);

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int iat2 = this->iwt2iat.at(iwt2);

				const int irc = ic*ParaO.nrow + ir;

				if(iat1==iat2 && iwt1==iwt2)
				{
					ftmp.at(iat1).at(dim) += dm_VU_dSm.at(dim).at(irc);
				}	
			}//end ic
		}//end ir

		for(int iat=0; iat<ucell.nat; iat++)
		{
			complex<double> val = ftmp.at(iat).at(dim);
			MPI_Allreduce(&val, &ftmp.at(iat).at(dim), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

			this->force_dftu.at(iat).at(dim) = ftmp.at(iat).at(dim).real();
		}
	}

	return;
}


void DFTU_RELAX::cal_stress_k(const vector<vector<complex<double>>>& VU)
{
	TITLE("DFTU_RELAX", "cal_stress_k");

	const char transN = 'N', transT = 'T', transC='C';
	const int  one_int = 1;
	//const double alpha = 1.0, beta = 0.0;
	const complex<double> alpha(1.0,0.0), beta(0.0,0.0);
  const complex<double> zero(0.0,0.0);
  
	int count = 0;
	for(int dim1=0; dim1<3; dim1++)
	{
		for(int dim2=dim1; dim2<3; dim2++)
		{
			vector<complex<double>> dm_VU_sover(ParaO.nloc, zero);

			for(int ik=0; ik<kv.nks; ik++)
			{
				const int spin = kv.isk[ik];
				
				// The first term
				vector<complex<double>> stress_tmp(ParaO.nloc);

				//Calculate mat_tmp=dm*VU
				vector<complex<double>> mat_tmp(ParaO.nloc);

				pzgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_k.at(ik).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);

				pzgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					this->soverlap_k[ik][count], &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					// dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
          dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
				}

				// The second term
				ZEROS(VECTOR_TO_PTR(stress_tmp), ParaO.nloc);

				pzgemm_(&transN, &transC,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					this->soverlap_k[ik][count], &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					// dm_VU_sover.at(irc) += 0.5*stress_tmp.at(irc);
          dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
				}

			}//end ik

			complex<double> stmp(0.0, 0.0);
			for(int ir=0; ir<ParaO.nrow; ir++)
			{
				const int iwt1 = ParaO.MatrixInfo.row_set[ir];

				for(int ic=0; ic<ParaO.ncol; ic++)
				{
					const int iwt2 = ParaO.MatrixInfo.col_set[ic];

					const int irc = ic*ParaO.nrow + ir;

					if(iwt1==iwt2) stmp += dm_VU_sover.at(irc);

				}//end ic

			}//end ir

			double val = stmp.real();
			MPI_Allreduce(&val, &stress_dftu.at(dim1).at(dim2), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			// complex<double> tmp;
			// MPI_Allreduce(&stmp, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
						
			count++;
		}//end dim2
		
	}//end dim1
	
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(i>j) this->stress_dftu.at(i).at(j) = stress_dftu.at(j).at(i);
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			this->stress_dftu.at(i).at(j) *=  ucell.lat0 / ucell.omega;
		}
	}

	return;
}


void DFTU_RELAX::cal_force_gamma(const vector<vector<double>> &VU)
{
	TITLE("DFTU_RELAX", "cal_force_gamma");

	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;
	
	vector<vector<double>> dm_VU_dSm;
	dm_VU_dSm.resize(3);
	for(int dim=0; dim<3; dim++)
	{
		dm_VU_dSm.at(dim).resize(ParaO.nloc, 0.0);
	}

	vector<double> force_tmp(ParaO.nloc, 0.0);
	//Calculate mat_tmp=dm*VU
	vector<double> mat_tmp(ParaO.nloc, 0.0);
	
	for(int ik=0; ik<kv.nks; ik++)
	{
		const int spin = kv.isk[ik];

		for(int dim=0; dim<3; dim++)
		{
			// The first part
			if(dim==0) //dim=1,2 are same as dim=0
			{
				ZEROS(VECTOR_TO_PTR(mat_tmp), ParaO.nloc);

				pdgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_gamma.at(spin).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);
			}

			if(dim==0)
			{
				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					LM.DSloc_x, &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==1)
			{
				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					LM.DSloc_y, &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==2)
			{
				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					LM.DSloc_z, &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			
			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) += force_tmp.at(irc);
			}

			// The second part
			ZEROS(VECTOR_TO_PTR(force_tmp), ParaO.nloc);

			if(dim==0)
			{
				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LM.DSloc_x, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==1)
			{
				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LM.DSloc_y, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}
			else if(dim==2)
			{
				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LM.DSloc_z, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(force_tmp), &one_int, &one_int, ParaO.desc);
			}

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				dm_VU_dSm.at(dim).at(irc) -= force_tmp.at(irc);
			}

		}// end dim

	}//end ik

	vector<vector<double>> ftmp;
	ftmp.resize(ucell.nat);
	for(int iat=0; iat<ucell.nat; iat++)
	{
		ftmp.at(iat).resize(3, 0.0);
	}

	for(int dim=0; dim<3; dim++)
	{
		for(int ir=0; ir<ParaO.nrow; ir++)
		{
			const int iwt1 = ParaO.MatrixInfo.row_set[ir];
			const int iat1 = this->iwt2iat.at(iwt1);

			for(int ic=0; ic<ParaO.ncol; ic++)
			{
				const int iwt2 = ParaO.MatrixInfo.col_set[ic];
				const int iat2 = this->iwt2iat.at(iwt2);

				const int irc = ic*ParaO.nrow + ir;

				if(iat1==iat2 && iwt1==iwt2)
				{
					ftmp.at(iat1).at(dim) += dm_VU_dSm.at(dim).at(irc);
				}	
			}//end ic
		}//end ir
	}//end dim


	for(int iat=0; iat<ucell.nat; iat++)
	{
		for(int dim=0; dim<3; dim++)
		{
			double tmp = 0.0;
			double val = ftmp.at(iat).at(dim);
			MPI_Allreduce(&val, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			this->force_dftu.at(iat).at(dim) = tmp;
		}

	}

	return;
}


void DFTU_RELAX::cal_stress_gamma(const vector<vector<double>> &VU)
{
	TITLE("DFTU_RELAX", "cal_stress_gamma");

	const char transN = 'N', transT = 'T';
	const int  one_int = 1;
	const double alpha = 1.0, beta = 0.0;
	
	int count = 0;
	for(int dim1=0; dim1<3; dim1++)
	{
		for(int dim2=dim1; dim2<3; dim2++)
		{
			vector<double> dm_VU_sover(ParaO.nloc, 0.0);

			for(int ik=0; ik<kv.nks; ik++)
			{
				const int spin = kv.isk[ik];
				
				// The first term
				vector<double> stress_tmp(ParaO.nloc, 0.0);

				//Calculate mat_tmp=dm*VU
				vector<double> mat_tmp(ParaO.nloc, 0.0);

				pdgemm_(&transT, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					LOC.wfc_dm_2d.dm_gamma.at(spin).c, &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(VU.at(spin)), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc);

				pdgemm_(&transN, &transN,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc, 
					this->soverlap_gamma[count], &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
				}

				// The second term
				ZEROS(VECTOR_TO_PTR(stress_tmp), ParaO.nloc);

				pdgemm_(&transN, &transT,
					&NLOCAL, &NLOCAL, &NLOCAL,
					&alpha, 
					this->soverlap_gamma[count], &one_int, &one_int, ParaO.desc, 
					VECTOR_TO_PTR(mat_tmp), &one_int, &one_int, ParaO.desc,
					&beta,
					VECTOR_TO_PTR(stress_tmp), &one_int, &one_int, ParaO.desc);

				for(int irc=0; irc<ParaO.nloc; irc++)
				{
					dm_VU_sover.at(irc) -= 0.5*stress_tmp.at(irc);
				}

			}//end ik

			double stmp = 0.0;
			for(int ir=0; ir<ParaO.nrow; ir++)
			{
				const int iwt1 = ParaO.MatrixInfo.row_set[ir];

				for(int ic=0; ic<ParaO.ncol; ic++)
				{
					const int iwt2 = ParaO.MatrixInfo.col_set[ic];

					const int irc = ic*ParaO.nrow + ir;

					if(iwt1==iwt2) stmp += dm_VU_sover.at(irc);

				}//end ic

			}//end ir

			MPI_Allreduce(&stmp, &stress_dftu.at(dim1).at(dim2), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
			
			count++;
		}//end dim2
		
	}//end dim1

	
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(i>j) this->stress_dftu.at(i).at(j) = stress_dftu.at(j).at(i);
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			this->stress_dftu.at(i).at(j) *=  ucell.lat0 / ucell.omega;
		}
	}

	return;
}


void DFTU_RELAX::folding_dSm_soverlap()
{
	TITLE("DFTU_RELAX", "folding_dSm_soverlap");

	int nnr = 0;

	if(GAMMA_ONLY_LOCAL)
	{
		if(STRESS)
		{
			for(int i=0; i<6; i++)
			{
				ZEROS(soverlap_gamma[i], ParaO.nloc);
			}			
		}
	}
	else
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			for(int dim=0; dim<3; dim++)
			{
				ZEROS(dSm_k[ik][dim], ParaO.nloc);
			}
		}

		if(STRESS)
		{
			for(int ik=0; ik<kv.nks; ik++)
			{
				for(int i=0; i<6; i++)
				{
					ZEROS(soverlap_k[ik][i], ParaO.nloc);
				}
			}
		}
	}
	

	  Vector3<double> tau1, tau2, dtau;
	  Vector3<double> dtau1, dtau2, tau0;
    for(int T1=0; T1<ucell.ntype; ++T1)
    {
		  Atom* atom1 = &ucell.atoms[T1];
      for(int I1=0; I1<atom1->na; ++I1)
      {
			  tau1 = atom1->tau[I1];
        const int start1 = ucell.itiaiw2iwt(T1,I1,0);    

        GridD.Find_atom(tau1, T1, I1);
        for(int ad=0; ad<GridD.getAdjacentNum()+1; ++ad)
        {
          const int T2 = GridD.getType(ad);
				  const int I2 = GridD.getNatom(ad);
          const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

				  Atom* atom2 = &ucell.atoms[T2];

				  tau2 = GridD.getAdjacentTau(ad);
				  dtau = tau2 - tau1;

				  double distance = dtau.norm() * ucell.lat0;
				  double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

          bool adj = false;
				  if(distance < rcut) adj = true;
				  else if(distance >= rcut)
				  {
				  	for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
				  	{
				  		const int T0 = GridD.getType(ad0); 
				  		const int I0 = GridD.getNatom(ad0); 
				  		const int iat0 = ucell.itia2iat(T0, I0);
				  		const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

				  		tau0 = GridD.getAdjacentTau(ad0);
				  		dtau1 = tau0 - tau1;
				  		dtau2 = tau0 - tau2;

				  		double distance1 = dtau1.norm() * ucell.lat0;
				  		double distance2 = dtau2.norm() * ucell.lat0;

				  		double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
				  		double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

				  		if( distance1 < rcut1 && distance2 < rcut2 )
				  		{
				  			adj = true;
				  			break;
				  		}
				  	}
				  }				

				  if(adj)
				  {
				  	for(int jj=0; jj<atom1->nw*NPOL; ++jj)
				  	{
              const int jj0 = jj/NPOL;

              const int iw1_all = start1 + jj0; 
              const int mu = ParaO.trace_loc_row[iw1_all];
					    if(mu<0)continue;

				  		const int L1 = atom1->iw2l[jj0];
				  		const int N1 = atom1->iw2n[jj0];
				  		const int m1 = atom1->iw2m[jj0];


				  		for(int kk=0; kk<atom2->nw*NPOL; ++kk)
				  		{
                const int kk0 = kk/NPOL;

                const int iw2_all = start2 + kk0;
						    const int nu = ParaO.trace_loc_col[iw2_all];
						    if(nu<0)continue;

				  			const int L2 = atom2->iw2l[kk0];
				  			const int N2 = atom2->iw2n[kk0];
				  			const int m2 = atom2->iw2m[kk0];
  
				  			// if ( !ParaO.in_this_processor(iw1_all,iw2_all) )
				  			// {
				  				// ++iw2_all;
				  				// continue;
				  			// }

				  			// int mu = ParaO.trace_loc_row[iw1_all];
				  			// int nu = ParaO.trace_loc_col[iw2_all];
				  			int irc = nu*ParaO.nrow + mu;
  
				  			if(GAMMA_ONLY_LOCAL)
							  {
							  	if(STRESS)
							  	{
							  		this->soverlap_gamma[0][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+0];
							  		this->soverlap_gamma[1][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+1];
							  		this->soverlap_gamma[2][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+2];
							  		this->soverlap_gamma[3][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+1];
							  		this->soverlap_gamma[4][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+2];
							  		this->soverlap_gamma[5][irc] += LM.DSloc_Rz[nnr]*LM.DH_r[nnr*3+2];
							  	}
							  }
				  			else
				  			{
				  				Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z); 
  
				  				for(int ik=0; ik<kv.nks; ik++)
				  				{
				  					const double arg = ( kv.kvec_d[ik] * dR ) * TWO_PI;
				  					const complex<double> kphase( cos(arg),  sin(arg) );

				  					this->dSm_k[ik][0][irc] += LM.DSloc_Rx[nnr]*kphase;
				  					this->dSm_k[ik][1][irc] += LM.DSloc_Ry[nnr]*kphase;
				  					this->dSm_k[ik][2][irc] += LM.DSloc_Rz[nnr]*kphase;

				  					if(STRESS)
				  					{		
				  						this->soverlap_k[ik][0][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+0]*kphase;
				  						this->soverlap_k[ik][1][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+1]*kphase;
				  						this->soverlap_k[ik][2][irc] += LM.DSloc_Rx[nnr]*LM.DH_r[nnr*3+2]*kphase;
				  						this->soverlap_k[ik][3][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+1]*kphase;
				  						this->soverlap_k[ik][4][irc] += LM.DSloc_Ry[nnr]*LM.DH_r[nnr*3+2]*kphase;
				  						this->soverlap_k[ik][5][irc] += LM.DSloc_Rz[nnr]*LM.DH_r[nnr*3+2]*kphase;																
				  					}
				  				}
				  			}
				  			++nnr;
				  		}// kk
				    }// jj
				  }// adj
				  // else if(distance>=rcut)
				  // {
				  	// int start1 = ucell.itiaiw2iwt( T1, I1, 0);
				  	// int start2 = ucell.itiaiw2iwt( T2, I2, 0);
				  	// bool is_adj = false;
				  	// for (int ad0=0; ad0<GridD.getAdjacentNum()+1; ++ad0)
				  	// {
				  	// 	const int T0 = GridD.getType(ad0);
				  		
				  	// 	tau0 = GridD.getAdjacentTau(ad0);
				  	// 	dtau1 = tau0 - tau1;
				  	// 	double distance1 = dtau1.norm() * ucell.lat0;
				  	// 	double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
				  	// 	dtau2 = tau0 - tau2;
				  	// 	double distance2 = dtau2.norm() * ucell.lat0;
				  	// 	double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();
				  	// 	if(distance1<rcut1 && distance2<rcut2)
				  	// 	{
				  	// 		is_adj = true;
				  	// 		break;
				  	// 	}
				  	// }//ad0
				  	// if( is_adj )
				  	// {
				  // 		for(int jj=0; jj<atom1->nw * NPOL; ++jj)
				  // 		{
				  // 			const int mu = ParaO.trace_loc_row[start1+jj];
				  // 			if(mu<0) continue; 

				  // 			for(int kk=0; kk<atom2->nw * NPOL; ++kk)
				  // 			{
				  // 				const int nu = ParaO.trace_loc_col[start2+kk];
				  // 				if(nu<0) continue;

				  // 				++nnr;
				  // 			}//kk
				  // 		}//jj
				  // 	// }
				  // }//distance
			  }// ad
		  }// I1
	  }// T1

	return;
}


void DFTU_RELAX::allocate_force_stress()
{
	TITLE("DFTU_RELAX","allocate_force_stress");

	if(GAMMA_ONLY_LOCAL)
	{
		if(STRESS)
		{
			// this->soverlap_gamma.resize(6);   //xx, xy, xz, yy, yz, zz
			// for(int i=0; i<6; i++)
			// {				
				// this->soverlap_gamma.at(i).resize(ParaO.nloc, 0.0);
			// }
			soverlap_gamma = new double* [6];
			for(int i=0; i<6; i++)
			{				
				soverlap_gamma[i] = new double [ParaO.nloc];
			}
		}
	}
	else
	{
		dSm_k = new complex<double>** [kv.nks];
		//this->dSm_k.resize(kv.nks);

		for(int ik=0; ik<kv.nks; ik++)
		{
			//this->dSm_k.at(ik).resize(3);
			dSm_k[ik] = new complex<double>* [3];
			for(int dim=0; dim<3; dim++)
			{
				//this->dSm_k.at(ik).at(dim).resize(ParaO.nloc, ZERO);
				dSm_k[ik][dim] = new complex<double>[ParaO.nloc];
			}
		}

		if(STRESS)
		{
			//this->soverlap_k.resize(kv.nks);
			soverlap_k = new complex<double>** [kv.nks];

			for(int ik=0; ik<kv.nks; ik++)
			{
				//this->soverlap_k.at(ik).resize(6);   //xx, xy, xz, yy, yz, zz
				soverlap_k[ik] = new complex<double>* [6];
				for(int i=0; i<6; i++)
				{				
					//this->soverlap_k.at(ik).at(i).resize(ParaO.nloc, ZERO);
					soverlap_k[ik][i] = new complex<double> [ParaO.nloc];
				}
			}
		}
	}
	
	return;
}


void DFTU_RELAX::erase_force_stress()
{
	TITLE("DFTU_RELAX","erase_force_stress");

	if(GAMMA_ONLY_LOCAL)
	{
		if(STRESS)
		{
			//this->soverlap_gamma.resize(6);   //xx, xy, xz, yy, yz, zz
			for(int i=0; i<6; i++)
			{				
				//this->soverlap_gamma.at(i).resize(1, 0.0);
				delete [] soverlap_gamma[i];
			}
			delete [] soverlap_gamma;
      soverlap_gamma=nullptr;
		}
	}
	else
	{
		//this->dSm_k.resize(kv.nks);

		for(int ik=0; ik<kv.nks; ik++)
		{
			//this->dSm_k.at(ik).resize(3);
			for(int dim=0; dim<3; dim++)
			{
				//this->dSm_k.at(ik).at(dim).resize(1, ZERO);
				delete [] dSm_k[ik][dim];
			}
			delete [] dSm_k[ik];
		}
		delete [] dSm_k;
    dSm_k = nullptr;

		if(STRESS)
		{
			//this->soverlap_k.resize(kv.nks);

			for(int ik=0; ik<kv.nks; ik++)
			{
				//this->soverlap_k.at(ik).resize(6);   //xx, xy, xz, yy, yz, zz
				for(int i=0; i<6; i++)
				{		
					//this->soverlap_k.at(ik).at(i).resize(1, ZERO);
					delete [] soverlap_k[ik][i];
				}
				delete [] soverlap_k[ik];
			}
			delete [] soverlap_k;
      soverlap_k = nullptr;
		}
	}
			
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
	TITLE("DFTU_RELAX","get_onebody_eff_pot");

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
