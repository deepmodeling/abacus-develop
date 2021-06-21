#include"stress_func.h"

//calculate the kinetic stress in PW base
void Stress_Func::stress_kin(matrix& sigma)
{
	double *kfac;
	double **gk;
	gk=new double* [3];
	double tbsp,gk2,arg;
	int ik,l,m,i,j,ibnd,is;
	int npw;
	double s_kin[3][3];
	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			s_kin[l][m]=0.0;
		}
	}
		
	tbsp=2.0/sqrt(PI);
		
	int npwx=0;
	int qtot = 0;
	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int ig=0; ig<kv.ngk[ik]; ig++)
		{
			qtot += kv.ngk[ik];
		}
		if(npwx<kv.ngk[ik])npwx=kv.ngk[ik];
	}
		
	kfac=new double[npwx];
	gk[0]= new double[npwx]; 
	gk[1]= new double[npwx];
	gk[2]= new double[npwx];
	for(i=0;i<npwx;i++)
	{
		kfac[i]=1;
	}
	double factor=TWO_PI/ucell.lat0;
	//    if(nks>1){
	//       iunigh.clear();
	//       iunigh.seekg(0,ios::beg);
	//    }//go back to the beginning of file

	for(ik=0;ik<kv.nks;ik++)
	{

		npw = kv.ngk[ik];
	//       if(kv.nks>1){
	//          iunigk>>igk;
	//          get_buffer(evc,nwordwfc,iunwfc,ik);
	//       }
		for(i=0;i<npw;i++)
		{
			gk[0][i] = pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, i), 0) * factor;
			gk[1][i] = pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, i), 1) * factor;
			gk[2][i] = pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, i), 2) * factor;

			//          if(qcutz>0){
			//             gk2=pow(gk[i].x,2)+pow(gk[i].y,2)+pow(gk[i].z,2);
			//             arg=pow((gk2-ecfixed)/q2sigma,2);
			//             kfac[i]=1+qcutz/q2sigma*tbsp*exp(-arg);
			//          }
		}

		//kinetic contribution

		for(l=0;l<3;l++)
		{
			for(m=0;m<l+1;m++)
			{
				for(ibnd=0;ibnd<NBANDS;ibnd++)
				{
					for(i=0;i<npw;i++)
					{
						if(0)
						{
							s_kin[l][m] +=
								wf.wg(ik,ibnd)*gk[l][i]*gk[m][i]*kfac[i]
								*(double((conj(wf.evc[ik](ibnd, i))
								*wf.evc[ik](ibnd, i)).real())+
								double((conj(wf.evc[ik](ibnd, i))*wf.evc[ik](ibnd, i+npwx)).real()));
						}
						else
						{
							s_kin[l][m] +=
								wf.wg(ik, ibnd)*gk[l][i]*gk[m][i]*kfac[i]
								*(double((conj(wf.evc[ik](ibnd, i))*wf.evc[ik](ibnd, i)).real()));
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
		
	for(l=0;l<3;l++)
	{
		for(m=0;m<l;m++)
		{
			s_kin[m][l]=s_kin[l][m];
		}
	}

	if(INPUT.gamma_only)
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				s_kin[l][m] *= 2.0*e2/ucell.omega;
			}
		}
	}
	else 
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				s_kin[l][m] *= e2/ucell.omega;
			}
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( s_kin[l][m] );
		}
	}


	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigma(l,m) = s_kin[l][m];
		}
	}
	//do symmetry
	if(Symmetry::symm_flag)
	{
		symm.stress_symmetry(sigma);
	}//end symmetry
	
	delete[] kfac;
	delete[] gk[0];
	delete[] gk[1];
	delete[] gk[2];
	delete[] gk;
		
	return;
}
