#include"stress_func.h"
#include "global.h"

//calculate the kinetic stress in PW base
void Stress_Func::stress_kin(ModuleBase::matrix& sigma)
{
	double *kfac;
	double **gk;
	gk=new double* [3];
	int npw;
	double s_kin[3][3];
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			s_kin[l][m]=0.0;
		}
	}
		
	int npwx=0;
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		if(npwx<GlobalC::kv.ngk[ik])npwx=GlobalC::kv.ngk[ik];
	}
		
	kfac=new double[npwx];
	gk[0]= new double[npwx]; 
	gk[1]= new double[npwx];
	gk[2]= new double[npwx];
	for(int i=0;i<npwx;i++)
	{
		kfac[i]=1;
	}
	double factor=ModuleBase::TWO_PI/GlobalC::ucell.lat0;

	for(int ik=0;ik<GlobalC::kv.nks;ik++)
	{
		npw = GlobalC::kv.ngk[ik];
		for(int i=0;i<npw;i++)
		{
			gk[0][i] = GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, i), 0) * factor;
			gk[1][i] = GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, i), 1) * factor;
			gk[2][i] = GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, i), 2) * factor;
		}

		//kinetic contribution

		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				for(int ibnd=0;ibnd<GlobalV::NBANDS;ibnd++)
				{
					for(int i=0;i<npw;i++)
					{
						if(0)
						{
							s_kin[l][m] +=
								GlobalC::wf.wg(ik,ibnd)*gk[l][i]*gk[m][i]*kfac[i]
								*(double((conj(GlobalC::wf.evc[ik](ibnd, i))
								*GlobalC::wf.evc[ik](ibnd, i)).real())+
								double((conj(GlobalC::wf.evc[ik](ibnd, i))*GlobalC::wf.evc[ik](ibnd, i+npwx)).real()));
						}
						else
						{
							s_kin[l][m] +=
								GlobalC::wf.wg(ik, ibnd)*gk[l][i]*gk[m][i]*kfac[i]
								*(double((conj(GlobalC::wf.evc[ik](ibnd, i))*GlobalC::wf.evc[ik](ibnd, i)).real()));
						}
					}
				}
			}
		}
		   
		//contribution from the nonlocal part
		   
		//stres_us(ik, gk, npw);
	}
		
	//add the US term from augmentation charge derivatives
		
	// addussstres(sigmanlc);
	
	//mp_cast
		
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l;m++)
		{
			s_kin[m][l]=s_kin[l][m];
		}
	}

	if(INPUT.gamma_only)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				s_kin[l][m] *= 2.0*ModuleBase::e2/GlobalC::ucell.omega;
			}
		}
	}
	else 
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				s_kin[l][m] *= ModuleBase::e2/GlobalC::ucell.omega;
			}
		}
	}

	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_all( s_kin[l][m] ); //qianrui fix a bug for kpar > 1
		}
	}


	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigma(l,m) = s_kin[l][m];
		}
	}
	//do symmetry
	if(ModuleSymmetry::Symmetry::symm_flag)
	{
		GlobalC::symm.stress_symmetry(sigma, GlobalC::ucell);
	}//end symmetry
	
	delete[] kfac;
	delete[] gk[0];
	delete[] gk[1];
	delete[] gk[2];
	delete[] gk;
		
	return;
}
