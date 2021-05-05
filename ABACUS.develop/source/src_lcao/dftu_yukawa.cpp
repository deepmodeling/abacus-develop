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

#include "dftu_yukawa.h"
#include "../src_global/constants.h"
#include "../src_pw/global.h"
#include "global_fp.h"
#include "../src_global/global_function.h"
#include "LOOP_ions.h"
#include "LCAO_matrix.h"


extern "C"
{
	void sphbsl_(int *n, double *r, double *A, double *val);
	void sphhnk_(int *n, double *r, double *A, double *val);
}

DFTU_Yukawa::DFTU_Yukawa(){}

DFTU_Yukawa::~DFTU_Yukawa(){}


void DFTU_Yukawa::cal_yukawa_lambda()
{
	TITLE("DFTU_Yukawa", "cal_yukawa_lambda");
	
	double sum_rho = 0.0;
	double sum_rho_lambda = 0.0;	
	for(int is=0; is<NSPIN; is++)
	{
		for(int ir=0; ir<pw.nrxx; ir++) 
		{
			double rho_ir = CHR.rho[is][ir];
			sum_rho += rho_ir;

			double lambda_ir = 2*pow(3*rho_ir/PI, (double)1.0/6.0);
			sum_rho_lambda += lambda_ir*rho_ir;
		}
	}

	double val1 = 0.0;
	double val2 = 0.0;
	MPI_Allreduce(&sum_rho, &val1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_rho_lambda, &val2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	this->lambda = val2/val1;

	//rescaling
	this->lambda /= 1.6;

	return;
}


void DFTU_Yukawa::cal_slater_Fk(const int L, const int T)
{
	TITLE("DFTU_Yukawa","cal_slater_Fk");
	
	if(Yukawa)
	{	 
		//this->lambda = INPUT.yukawa_lambda;
		
		for(int chi=0; chi<ucell.atoms[T].l_nchi[L]; chi++)
		{
		 //	if(chi!=0) continue;
			const int mesh = ORB.Phi[T].PhiLN(L,chi).getNr();
	
			for(int k=0; k<=L; k++)
			{			
				for(int ir0=1; ir0<mesh; ir0++)
				{
					double r0 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir0);
					const double rab0 = ORB.Phi[T].PhiLN(L,chi).getRab(ir0);
					const double R_L0 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir0);
	
					for(int ir1=1; ir1<mesh; ir1++) 
					{
						double bslval, hnkval;
						double r1 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir1);
						const double rab1 = ORB.Phi[T].PhiLN(L,chi).getRab(ir1);
						const double R_L1 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir1);		
						
						int l = 2*k;
						if(ir0<ir1)  //less than
						{
						 	sphbsl_(&l, &r0, &lambda, &bslval);
							sphhnk_(&l, &r1, &lambda, &hnkval);
						}
						else //greater than
						{
						 	sphbsl_(&l, &r1, &lambda, &bslval);
							sphhnk_(&l, &r0, &lambda, &hnkval);
						}					
						this->Fk.at(T).at(L).at(chi).at(k) -= (4*k+1)*lambda*pow(R_L0,2)*bslval*hnkval*pow(R_L1,2)*pow(r0,2)*pow(r1,2)*rab0*rab1;					
					}
				}
			 //	this->Fk.at(T).at(chi).at(k) /= pow(norm,2);
			}			
		}		
	}

	return;
}


void DFTU_Yukawa::cal_slater_UJ(const int istep, const int iter)
{
	TITLE("DFTU_Yukawa", "cal_slater_UJ");
	if(!Yukawa) return;

	this->cal_yukawa_lambda();
	
	for(int it=0; it<ucell.ntype; it++)
	{			
		const int NL = ucell.atoms[it].nwl + 1;
	
		for(int l=0; l<NL; l++)
		{
			int N = ucell.atoms[it].l_nchi[l];
			for(int n=0; n<N; n++)
			{
				ZEROS(VECTOR_TO_PTR(this->Fk.at(it).at(l).at(n)), l+1);
			}			
		}			 	
	}
 	
	for(int T=0; T<ucell.ntype; T++)
	{			
		const int NL = ucell.atoms[T].nwl + 1;

		for(int L=0; L<NL; L++)
		{
			const int N = ucell.atoms[T].l_nchi[L];

			if(L>=INPUT.orbital_corr[T] && INPUT.orbital_corr[T]!=-1)
			{
 if(L!=INPUT.orbital_corr[T]) continue;

				this->cal_slater_Fk(L, T);

				for(int n=0; n<N; n++)
				{
 if(n!=0) continue;

					switch(L)
					{
					case 1: //p electrons
    				    if(Yukawa)
						{
							this->U_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(0);
    				    	this->J_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(1)/5.0;
						}
						else
						{
						 //	if(n!=0) continue;

							this->Fk.at(T).at(L).at(n).at(0) = this->U[T];
    				    	this->Fk.at(T).at(L).at(n).at(1) = 5.0*this->J[T];
						}			
						break;

    				case 2: //d electrons
    				    if(Yukawa)
						{
							this->U_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(0);
    				    	this->J_Yukawa.at(T).at(L).at(n) = (this->Fk.at(T).at(L).at(n).at(1)+this->Fk.at(T).at(L).at(n).at(2))/14.0;
						}
						else
						{
						 //	if(n!=0) continue;

							this->Fk.at(T).at(L).at(n).at(0) = this->U[T];
							this->Fk.at(T).at(L).at(n).at(1) = 14.0*this->J[T]/(1.0+0.625);
							this->Fk.at(T).at(L).at(n).at(2) = 0.625*this->Fk.at(T).at(L).at(n).at(1);
						}
						break;

    				case 3: //f electrons
    				    if(Yukawa)
						{
							this->U_Yukawa.at(T).at(L).at(n) = this->Fk.at(T).at(L).at(n).at(0);
    				    	this->J_Yukawa.at(T).at(L).at(n) = (286.0*this->Fk.at(T).at(L).at(n).at(1) +
									195.0*this->Fk.at(T).at(L).at(n).at(2)+250.0*this->Fk.at(T).at(L).at(n).at(3))/6435.0;
						}
						else
						{
						 //	if(n!=0) continue;

							this->Fk.at(T).at(L).at(n).at(0) = this->U[T];
							this->Fk.at(T).at(L).at(n).at(1) = 6435.0*this->J[T]/(286.0+195.0*0.668+250.0*0.494);
							this->Fk.at(T).at(L).at(n).at(2) = 0.668*this->Fk.at(T).at(L).at(n).at(1);
							this->Fk.at(T).at(L).at(n).at(3) = 0.494*this->Fk.at(T).at(L).at(n).at(1);
						}
						break;
					}

					//Hartree to Rydeberg
					this->U_Yukawa.at(T).at(L).at(n) *= 2.0;
					this->J_Yukawa.at(T).at(L).at(n) *= 2.0;
				}//end n

			}//end if

		}//end L			 	
	}// end T
		
	return;
}

/*
void DFTU::cal_unscreened_slater_Fk(const int L, const int T)
{
	TITLE("DFTU","cal_slater_Fk");

	for(int chi=0; chi<ucell.atoms[T].l_nchi[L]; chi++)
	{
		const int mesh = ORB.Phi[T].PhiLN(L,chi).getNr();

		for(int k=0; k<=L; k++)
		{			
			for(int ir0=1; ir0<mesh; ir0++)
			{
				double r0 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir0);
				const double rab0 = ORB.Phi[T].PhiLN(L,chi).getRab(ir0);
				const double R_L0 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir0);

				for(int ir1=1; ir1<mesh; ir1++) 
				{
					double numerator, denominator;
					double r1 = ORB.Phi[T].PhiLN(L,chi).getRadial(ir1);
					const double rab1 = ORB.Phi[T].PhiLN(L,chi).getRab(ir1);
					const double R_L1 = ORB.Phi[T].PhiLN(L,chi).getPsi(ir1);		
					
					int l = 2*k;
					if(ir0<ir1)  //less than
					{
					 	numerator  = pow(r0, l);
						denominator = pow(r1, l+1);
					}
					else //greater than
					{
					 	numerator  = pow(r1, l);
						denominator = pow(r0, l+1);
					}					
					this->Fk.at(T).at(chi).at(k) += (numerator/denominator)*pow(R_L0,2)*pow(R_L1,2)*pow(r0,2)*pow(r1,2)*rab0*rab1;					
				}
			}
		}
		
		
	}

	this->cal_slater_UJ(T, L);

	return;
}

void DFTU::cal_slater_Vsc(const int T, const int L)
{
	TITLE("DFTU", "cal_slater_Vsc");

	for(int N=0; N<ucell.atoms[T].l_nchi[L]; N++)
	{
		if(!Yukawa && N!=0) continue;

		for(int m0=0; m0<2*L+1; m0++)
		{
			const int gindex0 = L*L + m0;

			for(int m1=0; m1<2*L+1; m1++)
			{				
				const int gindex1 = L*L + m1;

				for(int m2=0; m2<2*L+1; m2++)
				{
					const int gindex2 = L*L + m2;
					const int M0 = m0*(2*L+1) + m2;

					for(int m3=0; m3<2*L+1; m3++)
					{
						const int M1 = m1*(2*L+1) + m3;
						const int gindex3 = L*L + m3;

						for(int k=0; k<=2*L; k+=2)
						{
							const int l = (int)(k/2);
							for(int q=0; q<2*k+1; q++)
							{
								int gindex = k*k + q;

								double gaunt1 = UOT.get_Gaunt_coefficients(gindex0, gindex2, gindex);
								double gaunt2 = UOT.get_Gaunt_coefficients(gindex1, gindex3, gindex);

								this->Vsc.at(T).at(N)(M0, M1) += FOUR_PI*gaunt1*gaunt2*Fk.at(T).at(N).at(l)/(2.0*k+1.0);
							}
						}
					}
				}
			}
		}
	}

	return;
}
*/
