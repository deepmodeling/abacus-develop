#include "stress.h"
#include "global.h"
#include "potential.h"
#include "xc_functional.h"
#include "xc_gga_pw.h"
#include "efield.h"
#include "myfunc.h"
#include "symmetry.h"
// new
#include "H_Ewald_pw.h"
#include "H_Hartree_pw.h"
#include "H_XC_pw.h"

void Stress::cal_stress()
{
	TITLE("Stress","cal_stress");
	timer::tick("Stress","cal_stress",'E');    

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			sigmatot[i][j] = 0.0;
			sigmaxc[i][j] = 0.0;
			sigmahar[i][j] = 0.0;
			sigmakin[i][j] = 0.0;
			sigmaloc[i][j] = 0.0;
			sigmanlc[i][j] = 0.0;
			sigmaewa[i][j] = 0.0;
			sigmaxcc[i][j] = 0.0;
		}
	}
	//kinetic contribution
	stres_knl();

	if(Symmetry::symm_flag)
	{
		symm.stress_symmetry(sigmakin);
	}//end symmetry

	//hartree contribution
	stres_har();

    //ewald contribution
    stres_ewa();

    //xc contribution: add gradient corrections(non diagonal)
    for(int i=0;i<3;i++)
	{
       sigmaxc[i][i] = - (H_XC_pw::etxc - H_XC_pw::vtxc) / ucell.omega;
    }
    stres_gradcorr();

    //local contribution
    stres_loc();
    
    //nlcc
    stres_cc();
   
    //nonlocal
	stres_nl();

	if(Symmetry::symm_flag)
	{
		symm.stress_symmetry(sigmanlc);
	}//end symmetry

    for(int ipol=0;ipol<3;ipol++)
	{
        for(int jpol=0;jpol<3;jpol++)
		{
            sigmatot[ipol][jpol] = sigmakin[ipol][jpol] + sigmahar[ipol][jpol] + sigmanlc[ipol][jpol] 
                                 + sigmaxc[ipol][jpol] + sigmaxcc[ipol][jpol] + sigmaewa[ipol][jpol]
                                 + sigmaloc[ipol][jpol];
        }
    }
    
	if(Symmetry::symm_flag)
	{
		symm.stress_symmetry(sigmatot);
	}

/*    cout<<"total stress:"<<endl;
    for(int ipol=0;ipol<3;ipol++){
       for(int jpol=0;jpol<3;jpol++){
          cout << sigmatot(ipol,jpol)<<" ";
       }
       cout<< endl;
    }
  */ 
	bool ry = false;
	this->printstress_total(ry);

	if(TEST_STRESS) 
	{               
		ofs_running << "\n PARTS OF STRESS: " << endl;
		ofs_running << setiosflags(ios::showpos);
		ofs_running << setiosflags(ios::fixed) << setprecision(10) << endl;
		this->print_stress("KINETIC    STRESS",sigmakin,TEST_STRESS,ry);
		this->print_stress("LOCAL    STRESS",sigmaloc,TEST_STRESS,ry);
		this->print_stress("HARTREE    STRESS",sigmahar,TEST_STRESS,ry);
		this->print_stress("NON-LOCAL    STRESS",sigmanlc,TEST_STRESS,ry);
		this->print_stress("XC    STRESS",sigmaxc,TEST_STRESS,ry);
		this->print_stress("EWALD    STRESS",sigmaewa,TEST_STRESS,ry);
		this->print_stress("NLCC    STRESS",sigmaxcc,TEST_STRESS,ry);
		this->print_stress("TOTAL    STRESS",sigmatot,TEST_STRESS,ry);
	}
	return;
    
}


void Stress::print_stress(const string &name, double f[][3], const bool screen, bool ry)const
{
	ofs_running << " --------------------------- " << name << " ----------------------------" << endl;


	double fac = 1.0;

	if(!ry)
	{
		//     fac = Ry_to_eV / 0.529177;
	}

	cout << setprecision(8);
	cout << setiosflags(ios::showpos);

	if(screen)
	{
		cout << " ------------------- " << name << " --------------------" << endl;

	}

	const double output_acc = 1.0e-8;
	for (int i=0;i<3;i++)
	{
		ofs_running << setw(15)<< " ";
		if( abs(f[i][0]) >output_acc) ofs_running << setw(15) << f[i][0]*fac;
		else ofs_running << setw(15) << "0";
		if( abs(f[i][1]) >output_acc) ofs_running << setw(15) << f[i][1]*fac;
		else ofs_running << setw(15) << "0";
		if( abs(f[i][2]) >output_acc) ofs_running << setw(15) << f[i][2]*fac;
		else ofs_running << setw(15) << "0";
		ofs_running << endl;

		if(screen)
		{
			if( abs(f[i][0]) >output_acc) cout << setw(15) << f[i][0]*fac;
			else cout << setw(15) << "0";
			if( abs(f[i][1]) >output_acc) cout << setw(15) << f[i][1]*fac;
			else cout << setw(15) << "0";
			if( abs(f[i][2]) >output_acc) cout << setw(15) << f[i][2]*fac;
			else cout << setw(15) << "0";
			cout << endl;
		}


	}

	cout << resetiosflags(ios::showpos);

	return;
}

void Stress::printstress_total (bool ry)
{
// zhengdy update 2016-10-08
	double unit_transform = 1;

	if(!ry)
	{
		unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
	}
	//  cout.setf(ios::fixed);


	//ofs_running << setiosflags(ios::right);
	ofs_running << setprecision(8) << setiosflags(ios::showpos) << setiosflags(ios::fixed) << endl;
	NEW_PART("TOTAL-STRESS (KBAR)");//Ryd/(a.u.)^3

//        if(INPUT.stress_set == 1)
	if(TEST_STRESS)
	{
		stringstream ss;
		ss << global_out_dir << "STRESS.dat" ;
		ofstream ofs( ss.str().c_str() );
		if(!ofs)
		{
			cout << "open STRESS.dat error !" <<endl;
		}

		for(int i=0; i<3; i++)
		{
			ofs << "   " << sigmatot[i][0]*unit_transform 
			<< "   " << sigmatot[i][1]*unit_transform 
			<< "   " << sigmatot[i][2]*unit_transform << endl;
		}

		ofs.close();
	}

	if(TEST_STRESS)
	{
		cout << setiosflags(ios::fixed) << setprecision(6);
		cout << setiosflags(ios::showpos);
		cout << " ------------------- TOTAL      STRESS --------------------" << endl;
		cout << " " << setw(8) << "STRESS" << endl;
		ofs_running << " " << setw(12) << "STRESS" << endl;
	}
	for (int i=0; i<3; i++)
	{

		if(TEST_STRESS)
			cout << " " << setw(15) << sigmatot[i][0]*unit_transform << setw(15)
			<< sigmatot[i][1]*unit_transform << setw(15) << sigmatot[i][2]*unit_transform << endl;

		ofs_running << " " << setw(15) << sigmatot[i][0]*unit_transform << setw(15)
			<< sigmatot[i][1]*unit_transform << setw(15) << sigmatot[i][2]*unit_transform << endl;


	}
	ofs_running << setiosflags(ios::left);
	cout << resetiosflags(ios::showpos);

	return;
}


void Stress::stres_knl(void)
{
	TITLE("Stress","stres_nkl");

    double *kfac;
    double **gk;
    gk=new double* [3];
    double tbsp,gk2,arg;
    int ik,l,m,i,j,ibnd,is;
    int npw;
    
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
    for(i=0;i<npwx;i++){
       kfac[i]=1;
    }
    double factor=TWO_PI/ucell.lat0;
//    if(nks>1){
//       iunigh.clear();
//       iunigh.seekg(0,ios::beg);
//    }//go back to the beginning of file

    for(ik=0;ik<kv.nks;ik++){

       npw = kv.ngk[ik];
//       if(kv.nks>1){
//          iunigk>>igk;
//          get_buffer(evc,nwordwfc,iunwfc,ik);
//       }
       for(i=0;i<npw;i++){
          gk[0][i]=(kv.kvec_c[ik].x+pw.gcar[wf.igk(ik, i)].x)*factor;
          gk[1][i]=(kv.kvec_c[ik].y+pw.gcar[wf.igk(ik, i)].y)*factor;
          gk[2][i]=(kv.kvec_c[ik].z+pw.gcar[wf.igk(ik, i)].z)*factor;
          
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
						   sigmakin[l][m]=sigmakin[l][m]+
							   wf.wg(ik,ibnd)*gk[l][i]*gk[m][i]*kfac[i]
							   *(double((conj(wf.evc[ik](ibnd, i))
											   *wf.evc[ik](ibnd, i)).real())+
									   double((conj(wf.evc[ik](ibnd, i))*wf.evc[ik](ibnd, i+npwx)).real()));
					   }
					   else
					   {
						   sigmakin[l][m]=sigmakin[l][m]+
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
			sigmakin[m][l]=sigmakin[l][m];
		}
	}

	if(INPUT.gamma_only)
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				sigmakin[l][m]=2.0*e2/ucell.omega*sigmakin[l][m];
			}
		}
	}
	else 
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				sigmakin[l][m]=e2/ucell.omega*sigmakin[l][m];
			}
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigmakin[l][m] );
		}
	}

//    cout<<"stres_knl"<<endl;
//    for(l=0;l<3;l++){
//       for(m=0;m<3;m++){
  //        cout<<setprecision(3)<<setiosflags(ios::scientific)<<sigmakin[l][m]<<" ";
  //     }
    //   cout<<endl;
//    }
/*        cout<<"gk"<<endl;
    for(i=0;i<npw;i++){
        cout<<gk[0][i]<<" "<<gk[1][i]<<" "<<gk[2][i]<<endl;     
    } 
        cout<<"evc"<<endl;
    for(i=0;i<npw;i++){
       // cout<<wf.evc[ik](ibnd, i).real()<<endl;
    }*/
    //symmetrize stress
    //symmatrix();
    //symmatrix();

    delete[] kfac;
    delete[] gk[0];
    delete[] gk[1];
    delete[] gk[2];
    delete[] gk;
    
    return;
}



void Stress::stres_har(void)
{
	TITLE("Stress","stres_har");

	double shart,g2;
	const double eps=1e-8;
	int is,ig,l,m,nspin0;

	/*     for(l=0;l<3;l++){
		   for(m=0;m<3;m++){
		   sigmahar[l][m]=0.0;
		   }
		   }*/

	complex<double> *Porter = UFFT.porter;

	//  Hartree potential VH(r) from n(r)
	ZEROS( Porter, pw.nrxx );
	for(int is=0; is<NSPIN; is++)
	{
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			Porter[ir] += complex<double>( CHR.rho[is][ir], 0.0 );
		}
	}
	//=============================
	//  bring rho (aux) to G space
	//=============================
	pw.FFT_chg.FFT3D(Porter, -1);

	complex<double> *psic = new complex<double> [pw.nrxx];
	double *psic0 = new double[pw.nrxx];
	ZEROS( psic0, pw.nrxx);
	for(int is=0; is<NSPIN; is++)
	{
		daxpy (pw.nrxx, 1.0, CHR.rho[is],1,psic0,2 );
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			psic[ir] = complex<double>(psic0[ir], 0.0);
		}
	}

	pw.FFT_chg.FFT3D(psic, -1) ;


	double charge;
	if (pw.gstart == 1)
	{
		charge = ucell.omega * Porter[pw.ig2fftc[0]].real();
	}

	complex<double> *vh_g  = new complex<double>[pw.ngmc];
	ZEROS(vh_g, pw.ngmc);

	double g[3];
	//test
	//         int i=pw.gstart;
	//         cout<< "gstart " <<pw.gstart<<endl;
	double ehart=0;
	for (int ig = pw.gstart; ig<pw.ngmc; ig++)
	{
		const int j = pw.ig2fftc[ig];
		if(pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
		{
			const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);

			ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
			// vh_g[ig] = fac * Porter[j];

			shart= ( conj( Porter[j] ) * Porter[j] ).real()/(ucell.tpiba2 * pw.gg [ig]);

			g[0]=pw.gcar[ig].x;
			g[1]=pw.gcar[ig].y;
			g[2]=pw.gcar[ig].z;
			//test

			//cout<<"g "<<g[0]<<" "<<g[1]<<" "<<g[2]<<endl;
			//cout<<"gg "<<pw.gg[i]<<" "<<pw.gcar[i].norm2()<<endl;
			//cout<<"Porter "<<Porter[j]<<endl;
			//cout<<"tpiba2 "<<ucell.tpiba2<<endl;


			for(l=0;l<3;l++)
			{
				for(m=0;m<l+1;m++)
				{
					sigmahar[l][m]=sigmahar[l][m]+shart *2*g[l]*g[m]/pw.gg[ig];
				}
			}
		}

	}
	Parallel_Reduce::reduce_double_pool( ehart );

	ehart *= 0.5 * ucell.omega;

	//cout<<"ehart "<<ehart<<" en.ehart "<< en.ehart<<endl;
	for(l=0;l<3;l++)
	{
		for(m=0;m<l+1;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigmahar[l][m] );
		}
	}

	//        Parallel_Reduce::reduce_double_pool( ehart );
	//        ehart *= 0.5 * ucell.omega;

	//psic(:)=(0.0,0.0)
	/*        nspin0=NSPIN;
			  if(NSPIN==4)nspin0=1;
			  for(is=0;is<nspin0;is++){
			  daxpy();
			  }*/

	if(INPUT.gamma_only)
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				sigmahar[l][m] *= e2 * FOUR_PI;
			}
		}
	}
	else
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				sigmahar[l][m] *= 0.5 * e2 * FOUR_PI;
			}
		}
	}

	for(l=0;l<3;l++)
	{
		sigmahar[l][l] -= H_Hartree_pw::hartree_energy /ucell.omega;
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<l;m++)
		{
			sigmahar[m][l]=sigmahar[l][m];
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigmahar[l][m]*=-1;
		}
	}

	//     cout<<"stres_har"<<endl;
	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			//         cout<<setprecision(3)<<setiosflags(ios::scientific)<<sigmahar[l][m]<<" ";
		}
		//    cout<<endl;
	}

	delete[] vh_g;
	delete[] psic;
	delete[] psic0;
	return;
}

void Stress::stres_loc(void)
{
	TITLE("Stress","stres_loc");

	double *dvloc;
	double evloc,fact;
	int ng,nt,l,m,is;

	dvloc = new double[pw.ngmc];

	/*    for(l=0;l<3;l++){
		  for(m=0;m<3;m++){
		  sigmaloc[l][m]=0;
		  }
		  }*/


	complex<double> *Porter = UFFT.porter;

	ZEROS( Porter, pw.nrxx );
	for(int is=0; is<NSPIN; is++)
	{
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			Porter[ir] += complex<double>(CHR.rho[is][ir], 0.0 );
		}
	}
	pw.FFT_chg.FFT3D(Porter, -1);

	/*        complex<double> *psic = new complex<double> [pw.nrxx];
			  double *psic0 = new double[pw.nrxx];
			  ZEROS( psic0, pw.nrxx);
			  for(int is=0; is<NSPIN; is++)
			  {
			  daxpy (pw.nrxx, 1.0, CHR.rho[is],1,psic0,2 );
			  for (int ir=0; ir<pw.nrxx; ir++)
			  {
			  psic[ir] = complex<double>(psic0[ir], 0.0);
			  }
			  }

			  pw.FFT_chg.FFT3D(psic, -1) ;
	 */  

	if(INPUT.gamma_only==1) fact=2.0;
	else fact=1.0;

	evloc=0.0;
	double g[3]={0,0,0};

	complex<double> *vg = new complex<double>[pw.ngmc];
	ZEROS( vg, pw.ngmc );
	for (int it=0; it<ucell.ntype; it++)
	{
		if (pw.gstart==1) 
		{
			evloc += ppcell.vloc(it, pw.ig2ngg[0]) * (pw.strucFac(it,0) * conj(Porter[pw.ig2fftc[0]])).real();
		}
		for (int ig=pw.gstart; ig<pw.ngmc; ig++)
		{
			const int j = pw.ig2fftc[ig];
			evloc += ppcell.vloc(it, pw.ig2ngg[ig]) * (pw.strucFac(it,ig) * conj(Porter[j]) * fact).real();
		}
	}
	for( nt = 0;nt< ucell.ntype; nt++)
	{
		const Atom* atom = &ucell.atoms[nt];
		//mark by zhengdy for check
		// if ( ppcell.vloc == NULL ){
		if(0)
		{
			// special case: pseudopotential is coulomb 1/r potential
			dvloc_coul (atom->zv, dvloc);
		}
		else{
			// normal case: dvloc contains dV_loc(G)/dG
			dvloc_of_g ( atom->msh, atom->rab, atom->r,
					atom->vloc_at, atom->zv, dvloc);
		}

		for( ng = 0;ng< pw.ngmc;ng++)
		{
			const int j = pw.ig2fftc[ng];
			g[0]=pw.gcar[ng].x;
			g[1]=pw.gcar[ng].y;
			g[2]=pw.gcar[ng].z;
			for (l = 0;l< 3;l++)
			{
				for (m = 0; m<l+1;m++)
				{
					sigmaloc[l][m] = sigmaloc[l][m] +  ( conj( Porter[j] ) 
							* pw.strucFac (nt, ng) ).real() * 2.0 * dvloc [ pw.ig2ngg[ng] ] 
						* ucell.tpiba2 * g[l] * g[m] * fact;
				}
			}
		}
	}

	for( l = 0;l< 3;l++)
	{
		sigmaloc [l][l] = sigmaloc [l][l] + evloc;
		for (m = 0; m<l+1; m++)
		{
			sigmaloc [m][l] = sigmaloc [l][m];
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigmaloc[l][m] );
		}
	}

	//mp_sum(  sigmaloc, intra_bgrp_comm )
	delete[] dvloc;
	delete[] vg; 
	return;
}
    
void Stress::dvloc_of_g (const int& msh,
                         const double* rab,
                         const double* r,
                         const double* vloc_at,
                         const double& zp,
                         double*  dvloc){
  //----------------------------------------------------------------------
  //
  // dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  //

  //
  //    first the dummy variables
  //
  //int ngl, mesh, msh;
      // the number of shell of G vectors
      // max number of mesh points
      // number of mesh points for radial integration

  //double zp, rab [mesh], r [mesh], vloc_at [mesh], tpiba2, omega, gl [ngl];
      // valence pseudocharge
      // the derivative of the radial grid
      // the radial grid
      // the pseudo on the radial grid
      // 2 pi / alat
      // the volume of the unit cell
      // the moduli of g vectors for each s
  //
  //double  dvloc[ngl];
  // the fourier transform dVloc/dG
  //
  double vlcp=0;
  double  *aux, *aux1;

  int i, igl, igl0;
  // counter on erf functions or gaussians
  // counter on g shells vectors
  // first shell with g != 0

  aux = new double[msh];
  aux1 = new double[msh];
  ZEROS(aux, msh);
  ZEROS(aux1, msh);

  // the  G=0 component is not computed
  if (pw.ggs[0] < 1.0e-8){
     dvloc[0] = 0.0;
     igl0 = 1;
  }
  else{
     igl0 = 0;
  }

  // Pseudopotentials in numerical form (Vloc contains the local part)
  // In order to perform the Fourier transform, a term erf(r)/r is
  // subtracted in real space and added again in G space

  //
  //   This is the part of the integrand function
  //   indipendent of |G| in real space
  //
  for( i = 0;i< msh; i++){
     aux1[i] = r [i] * vloc_at [i] + zp * e2 * erf(r[i]);
  }
  for( igl = igl0;igl< pw.nggm;igl++){
     double gx = sqrt (pw.ggs [igl] * ucell.tpiba2);
     double gx2 = pw.ggs [igl] * ucell.tpiba2;
     //
     //    and here we perform the integral, after multiplying for the |G|
     //    dependent  part
     //
     // DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
     for( i = 1;i< msh;i++){
        aux [i] = aux1 [i] * (r [i] * cos (gx * r [i] ) / gx - sin (gx * r [i] ) / pow(gx,2));
     }
    // simpson (msh, aux, rab, vlcp);
     Mathzone::Simpson_Integral(msh, aux, rab, vlcp );
     // DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
     vlcp *= FOUR_PI / ucell.omega / 2.0 / gx;
     // subtract the long-range term
     double g2a = gx2 / 4.0;
     vlcp += FOUR_PI / ucell.omega * zp * e2 * exp ( - g2a) * (g2a + 1) / pow(gx2 , 2);
     dvloc [igl] = vlcp;
  }
  delete[] aux1;
  delete[] aux;

  return;
}

void Stress::dvloc_coul (const double& zp,
                         double* dvloc){
  //----------------------------------------------------------------------
  //
  //    Fourier transform of the Coulomb potential - For all-electron
  //    calculations, in specific cases only, for testing purposes
  //

  //integer, intent(in) :: ngl
        // the number of shell of G vectors
  //real(DP), intent(in) :: zp, tpiba2, omega, gl (ngl)
        // valence pseudocharge
        // 2 pi / alat
        // the volume of the unit cell
        // the moduli of g vectors for each s
  //real(DP), intent(out) :: dvloc (ngl)
        // fourier transform: dvloc = D Vloc (g^2) / D g^2 = 4pi e^2/omegai /G^4

  int  igl0;
  // first shell with g != 0

  // the  G=0 component is 0
  if (pw.ggs[0] < 1.0e-8){
     dvloc[0] = 0.0;
     igl0 = 1;
  }
  else{
     igl0 = 0;
  }
  for(int i=igl0;i<pw.nggm;i++)
  dvloc[i] = FOUR_PI * zp * e2 / ucell.omega / pow(( ucell.tpiba2 * pw.ggs[i] ),2);

  return;
}

       

void Stress::stres_nl(){
	TITLE("Stress","stres_nl");
	timer::tick("Stress","stres_nl");


	const int nkb = ppcell.nkb;
	if(nkb == 0) return; 

	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	ComplexMatrix dbecp( nkb, NBANDS);
	ComplexMatrix becp( nkb, NBANDS);

	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ComplexMatrix vkb1( nkb, wf.npwx );
	ComplexMatrix vkb0[3];
	for(int i=0;i<3;i++){
		vkb0[i].create(nkb, wf.npwx);
	}
	ComplexMatrix vkb2( nkb, wf.npwx );
	for (int ik = 0;ik < kv.nks;ik++)
	{
		for(int i=0;i<3;i++){
			vkb0[i].zero_out();
		}
		vkb2.zero_out();      

		if (NSPIN==2) CURRENT_SPIN = kv.isk[ik];
		wf.npw = kv.ngk[ik];
		// generate vkb
		if (ppcell.nkb > 0)
		{
			ppcell.getvnl(ik);
		}

		// get becp according to wave functions and vkb
		// important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
		becp.zero_out();
		for (int ib=0; ib<NBANDS; ib++)
		{
			for (int i=0;i<nkb;i++)
			{
				for (int ig=0; ig<wf.npw; ig++)
				{
					becp(i,ib) += wf.evc[ik](ib,ig)* conj( ppcell.vkb(i,ig) );
				}
			}
		}
		Parallel_Reduce::reduce_complex_double_pool( becp.c, becp.size);



		for(int i=0;i<3;i++) {
			get_dvnl1(vkb0[i],ik,i);
		}

		get_dvnl2(vkb2,ik);

		Vector3<double> qvec;
		double qvec0[3];

		for (int ipol = 0; ipol<3; ipol++)
		{
			for(int jpol = 0; jpol < ipol+1; jpol++)
			{
				dbecp.zero_out();    
				vkb1.zero_out();
				for (int i = 0;i < nkb;i++)
				{

					for (int ig=0; ig<wf.npw; ig++)  
					{

						qvec = wf.get_1qvec_cartesian(ik,ig) ;
						qvec0[0] = qvec.x;
						qvec0[1] = qvec.y;
						qvec0[2] = qvec.z;


						vkb1(i, ig) += 0.5 * qvec0[ipol] * vkb0[jpol](i,ig)
							+ 0.5 * qvec0[jpol] * vkb0[ipol](i,ig) ;
					}//end ig


					/*     if (kpol==0)
						   {
						   for (int ig=0; ig<wf.npw; ig++)
						   {
						   qvec = wf.get_1qvec_cartesian(ik,ig);
						   vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * qvec.x;
						   }
						   }
						   if (kpol==1)
						   {
						   for (int ig=0; ig<wf.npw; ig++)
						   {
						   qvec = wf.get_1qvec_cartesian(ik,ig);
						   vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * qvec.y;
						   }
						   }
						   if (kpol==2)
						   {
						   for (int ig=0; ig<wf.npw; ig++)
						   {
						   qvec = wf.get_1qvec_cartesian(ik,ig);
						   vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * qvec.z;
						   }
						   }*/

				}//end nkb


				for (int ib=0; ib<NBANDS; ib++)
				{
					for (int i=0; i<nkb; i++)
					{
						for (int ig=0; ig<wf.npw; ig++)
						{
							//first term
							dbecp(i,ib) = dbecp(i,ib) - 2.0 * wf.evc[ik](ib,ig) * conj( vkb1(i,ig) ) ;
							//second termi
							if(ipol == jpol)
								dbecp(i,ib) += -1.0 * wf.evc[ik](ib,ig)* conj( ppcell.vkb(i,ig) );
							//third term
							qvec = wf.get_1qvec_cartesian(ik,ig);
							qvec0[0] = qvec.x;
							qvec0[1] = qvec.y;
							qvec0[2] = qvec.z;
							double qm1; 
							if(qvec.norm() > 1e-8) qm1 = 1.0 / qvec.norm();
							else qm1 = 0;
							dbecp(i,ib) +=  -2.0 * wf.evc[ik](ib,ig) * conj(vkb2(i,ig)) * qvec0[ipol] * qvec0[jpol] * qm1 * ucell.tpiba;
						}//end ig
					}//end i
				}//end ib


				//              don't need to reduce here, keep dbecp different in each processor,
				//              and at last sum up all the forces.
				//              Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

				//              double *cf = new double[ucell.nat*3];
				//              ZEROS(cf, ucell.nat);
				for (int ib=0; ib<NBANDS; ib++)
				{
					double fac = wf.wg(ik, ib) * 1.0;
					int iat = 0;
					int sum = 0;
					for (int it=0; it<ucell.ntype; it++)
					{
						const int Nprojs = ucell.atoms[it].nh;
						for (int ia=0; ia<ucell.atoms[it].na; ia++)
						{
							for (int ip=0; ip<Nprojs; ip++)
							{
								double ps = ppcell.deeq(CURRENT_SPIN, iat, ip, ip) ;
								const int inkb = sum + ip;
								//out<<"\n ps = "<<ps;


								const double dbb = ( conj( dbecp( inkb, ib) ) * becp( inkb, ib) ).real();
								sigmanlc[ipol][ jpol] -= ps * fac * dbb;

							}//end ip
							++iat;        
							sum+=Nprojs;
						}//ia
					} //end it
				} //end band

			}//end jpol
		}//end ipol
	}// end ik

	// sum up forcenl from all processors
	for(int l=0;l<3;l++){
		for(int m=0;m<3;m++){
			if(m>l) sigmanlc[l][m] = sigmanlc[m][l];
			Parallel_Reduce::reduce_double_pool( sigmanlc[l][m] );
		}
	}

	//        Parallel_Reduce::reduce_double_all(sigmanl.c, sigmanl.nr * sigmanl.nc);

	//        cout<<"sigmanl:"<<endl;
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigmanlc[ipol][jpol] *= 1.0 / ucell.omega;
			//                cout << sigmanl(ipol,jpol)<<" ";
		}
		//          cout<<endl;
	}

	//  this->print(ofs_running, "nonlocal stress", stresnl);
	timer::tick("Stress","stres_nl");
	return;
}
 
void Stress::get_dvnl1(
		ComplexMatrix &vkb,
		const int ik,
		const int ipol)
{
	if(test_pp) TITLE("Stress","get_dvnl1");
	//        timer::tick("Stress","get_dvnl1");

	const int lmaxkb = ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = kv.ngk[ik];
	const int nhm = ppcell.nhm;
	int ig, ia, nb, ih;
	matrix vkb1(nhm, npw);
	vkb1.zero_out();
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	matrix dylm(x1, npw);
	Vector3<double> *gk = new Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}

	dylmr2(x1, npw, gk, dylm, ipol);

	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		if(test_pp>1) OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = ucell.atoms[it].nbeta;
		const int nh = ucell.atoms[it].nh;

		if(test_pp>1) OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(test_pp>1) OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * ucell.tpiba;

				//                              cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
				//                              cout << "\n gk.norm = " << gnorm;

				vq [ig] = Mathzone::Polynomial_Interpolation(
						ppcell.tab, it, nb, NQX, DQ, gnorm );

			} // enddo

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = dylm(lm, ig) * vq [ig];
					}

				}

			} // end ih

		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<ucell.atoms[it].na; ia++)
		{
			complex<double> *sk = wf.get_sk(ik, it, ia);
			for (ih = 0;ih < nh;ih++)
			{
				complex<double> pref = pow( NEG_IMAG_UNIT, ppcell.nhtol(it, ih));      //?
				for (ig = 0;ig < npw;ig++)
				{
					vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	} // enddo
	delete [] gk;
	delete [] vq;
	//timer::tick("Stress","get_dvnl1");
	return;
}//end get_dvnl1


void Stress::get_dvnl2(
	ComplexMatrix &vkb,
	const int ik)
{
	if(test_pp) TITLE("Stress","get_dvnl2");
	timer::tick("Stress","get_dvnl2");

	const int lmaxkb = ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = kv.ngk[ik];
	const int nhm = ppcell.nhm;
	int ig, ia, nb, ih;
	matrix vkb1(nhm, npw);
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	matrix ylm(x1, npw);
	Vector3<double> *gk = new Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}
	Mathzone::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		if(test_pp>1) OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = ucell.atoms[it].nbeta;
		const int nh = ucell.atoms[it].nh;

		if(test_pp>1) OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(test_pp>1) OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * ucell.tpiba;

				//                              cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
				//                              cout << "\n gk.norm = " << gnorm;
				vq [ig] = Polynomial_Interpolation_nl(
						ppcell.tab, it, nb, DQ, gnorm );

			} // enddo

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<ucell.atoms[it].na; ia++)
		{
			complex<double> *sk = wf.get_sk(ik, it, ia);
			for (ih = 0;ih < nh;ih++)
			{
				complex<double> pref = pow( NEG_IMAG_UNIT, ppcell.nhtol(it, ih));      //?
				for (ig = 0;ig < npw;ig++)
				{
					vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
			delete [] sk;

		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;
	timer::tick("Stress","get_dvnl2");

	return;
}


double Stress::Polynomial_Interpolation_nl
(
    const realArray &table,
    const int &dim1,
    const int &dim2,
    const double &table_interval,
    const double &x                             // input value
)
{
//      timer::tick("Mathzone","Poly_Interpo_2");
    assert(table_interval>0.0);
    const double position = x  / table_interval;
    const int iq = static_cast<int>(position);

    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y=
      ( table(dim1, dim2, iq)   * (-x2*x3-x1*x3-x1*x2) / 6.0 +
        table(dim1, dim2, iq+1) * (+x2*x3-x0*x3-x0*x2) / 2.0 -
        table(dim1, dim2, iq+2) * (+x1*x3-x0*x3-x0*x1) / 2.0 +
        table(dim1, dim2, iq+3) * (+x1*x2-x0*x2-x0*x1) / 6.0 )/table_interval ;

//      timer::tick("Mathzone","Poly_Interpo_2");
    return y;
}

void Stress::dylmr2 (
         const int nylm,
         const int ngy,
         Vector3<double> *gk,
         matrix &dylm,
         const int ipol){
  //-----------------------------------------------------------------------
  //
  //     compute \partial Y_lm(G) \over \partial (G)_ipol
  //     using simple numerical derivation (SdG)
  //     The spherical harmonics are calculated in ylmr2
  //
  //int nylm, ngy, ipol;
  // number of spherical harmonics
  // the number of g vectors to compute
  // desired polarization
  //double g (3, ngy), gg (ngy), dylm (ngy, nylm)
  // the coordinates of g vectors
  // the moduli of g vectors
  // the spherical harmonics derivatives
  //
  int ig, lm;
  // counter on g vectors
  // counter on l,m component

  const double delta = 1e-6;
  double *dg, *dgi;

  matrix ylmaux;
  // dg is the finite increment for numerical derivation:
  // dg = delta |G| = delta * sqrt(gg)
  // dgi= 1 /(delta * sqrt(gg))
  // gx = g +/- dg


  Vector3<double> *gx = new Vector3<double> [ngy];
 

  dg = new double [ngy];
  dgi = new double [ngy];

  ylmaux.create (nylm, ngy);

  dylm.zero_out();
  ylmaux.zero_out();

  for( ig = 0;ig< ngy;ig++){
     gx[ig] = gk[ig];
  }
//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
  for( ig = 0;ig< ngy;ig++){
     dg [ig] = delta * gx[ig].norm() ;
     if (gx[ig].norm2() > 1e-9) {
        dgi [ig] = 1.0 / dg [ig];
     }
     else{
        dgi [ig] = 0.0;
     }
  }
//$OMP END PARALLEL DO

//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
  for( ig = 0;ig< ngy;ig++){
     if(ipol==0)
         gx [ig].x = gk[ ig].x + dg [ig];
     else if(ipol==1)
         gx [ig].y = gk [ ig].y + dg [ig];
     else if(ipol==2)
         gx [ig].z = gk [ ig].z + dg [ig];
  }
//$OMP END PARALLEL DO

  Mathzone::Ylm_Real(nylm, ngy, gx, dylm);
//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
  for(ig = 0;ig< ngy;ig++){
     if(ipol==0)
         gx [ig].x = gk [ ig].x - dg [ig];
     else if(ipol==1)
         gx [ig].y = gk [ ig].y - dg [ig];
     else if(ipol==2)
         gx [ig].z = gk [ ig].z - dg [ig];
  }
//$OMP END PARALLEL DO

  Mathzone::Ylm_Real(nylm, ngy, gx, ylmaux);


//  zaxpy ( - 1.0, ylmaux, 1, dylm, 1);
    for( lm = 0;lm< nylm;lm++){
        for(ig = 0;ig< ngy;ig++){
            dylm (lm,ig) = dylm(lm,ig) - ylmaux(lm,ig);
        }
    }


  for( lm = 0;lm< nylm;lm++){
//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
     for(ig = 0;ig< ngy;ig++){
        dylm (lm,ig) = dylm(lm,ig) * 0.5 * dgi [ig];
     }
//$OMP END PARALLEL DO
  }
  delete[] gx;
  delete[] dg;
  delete[] dgi;

  return;
}





//main part
void Stress::stres_ewa(){
    int i,j,l,m;
    double g[3];
    //the erfc function

    //tpiba2=pow(tpi/alat,2);
/*    for(i=0;i<3;i++){
       for(j=0;j<3;j++){
           sigmaewa[i][j]=0;
       }
    }*/
    double charge=0;
    for(int it=0; it < ucell.ntype; it++){
         for(int i=0; i<ucell.atoms[it].na; i++){
            charge = charge + ucell.atoms[it].zv;
         }
    }
    
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    double alpha=2.9;
    double upperbound;
    do
	{
       alpha-=0.1;
       if(alpha==0.0)
          WARNING_QUIT("stres_ew", "optimal alpha not found");
       upperbound =e2 * pow(charge,2) * sqrt( 2 * alpha / (TWO_PI)) * erfc(sqrt(ucell.tpiba2 * pw.ggchg / 4.0 / alpha));
    }
    while(upperbound>1e-7);
    
    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    double sdewald; 
    if(pw.gstart == 1)
	{
       sdewald = (TWO_PI) * e2 / 4.0 / alpha * pow(charge/ucell.omega,2);
    }
    else {
       sdewald = 0.0;
    }
    //sdewald is the diagonal term 

    double fact;
    if (INPUT.gamma_only) fact=2.0;
    else fact=1.0;
    
    double g2,g2a;
    double arg;
    complex<double> rhostar;
    double sewald;
    for(int ng=pw.gstart;ng<pw.ngmc;ng++){
        if(pw.gg[ng] >= 1.0e-12) //LiuXh 20180410
        {
            g2 = pw.gg[ng]* ucell.tpiba2;
            g2a = g2 /4.0/alpha;
            rhostar=(0.0,0.0);
            for(int it=0; it < ucell.ntype; it++){
                     for(int i=0; i<ucell.atoms[it].na; i++){
                            arg = (pw.gcar[ng].x * ucell.atoms[it].tau[i].x + pw.gcar[ng].y * ucell.atoms[it].tau[i].y + pw.gcar[ng].z * ucell.atoms[it].tau[i].z) * (TWO_PI);
                            rhostar = rhostar + complex<double>(ucell.atoms[it].zv * cos(arg),ucell.atoms[it].zv * sin(arg));
                     }
            }
            rhostar /= ucell.omega;
            sewald = fact* (TWO_PI) * e2 * exp(-g2a) / g2 * pow(abs(rhostar),2);
            sdewald = sdewald - sewald;
            g[0]=pw.gcar[ng].x;
            g[1]=pw.gcar[ng].y;
            g[2]=pw.gcar[ng].z;
            for(l=0;l<3;l++){
               for(m=0;m<l+1;m++){
                  sigmaewa[l][m]+=sewald * ucell.tpiba2 * 2.0 * g[l] * g[m] / g2 * (g2a+1);
               }
            }
        } //LiuXh 20180410
    }

    for(l=0;l<3;l++){
       sigmaewa[l][l]+=sdewald;
    }
    

    //R-space sum here (only for the processor that contains G=0) 
    int mxr = 50;
    int *irr;
    Vector3<double> *r;
    double *r2;
    r  = new Vector3<double>[mxr];
    r2 = new double[mxr];
    irr = new int[mxr];
    double rr;
    Vector3<double> d_tau;
    double r0[3];
    int rmax , nrm=0;
    double fac;
    if(pw.gstart==1){
       rmax = 4.0/sqrt(alpha)/ucell.lat0;
       //with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2*10^-1)
       for(int it=0; it < ucell.ntype; it++){
         for(int i=0; i<ucell.atoms[it].na; i++){
           for(int jt=0; jt < ucell.ntype; jt++){
             for(int j=0; j<ucell.atoms[jt].na; j++){
               //calculate tau[na]-tau[nb]
               d_tau = ucell.atoms[it].tau[i] - ucell.atoms[jt].tau[j];
               //generates nearest-neighbors shells 
               H_Ewald_pw::rgen(d_tau, rmax, irr, ucell.latvec, ucell.G, r, r2, nrm);
               for(int nr=0 ; nr<nrm ; nr++){
                rr=sqrt(r2[nr]) * ucell.lat0;
                fac = -e2/2.0/ucell.omega*pow(ucell.lat0,2)*ucell.atoms[it].zv * ucell.atoms[it].zv / pow(rr,3) * (erfc(sqrt(alpha)*rr)+rr * sqrt(8 * alpha / (TWO_PI)) * exp(-alpha * pow(rr,2)));
                for(l=0; l<3; l++){
                   for(m=0; m<l+1; m++){
                      r0[0] = r[nr].x;
                      r0[1] = r[nr].y;
                      r0[2] = r[nr].z;
                      sigmaewa[l][m] += fac * r0[l] * r0[m];
                   }//end m
                }//end l
             }//end nr
          }//end j
        }//end jt
       }//end i
      }//end it
    }//end if


    for(l=0;l<3;l++){
       for(m=0;m<l+1;m++){
          sigmaewa[m][l]=sigmaewa[l][m];
       }
    }
    for(l=0;l<3;l++){
       for(m=0;m<3;m++){
          sigmaewa[l][m]=-sigmaewa[l][m];
	  Parallel_Reduce::reduce_double_pool( sigmaewa[l][m] ); 
       }
    }
    
//    cout<<"stres_ewa: "<<endl;
//    for(l=0;l<3;l++){
//        for(m=0;m<3;m++){
//           cout<<setprecision(3)<<setiosflags(ios::scientific)<<sigmaewa[l][m]<<" ";
//        }
//        cout<<endl;
//    }
     
    delete [] r;
    delete [] r2;
    delete [] irr; 

    //mp_sum
    return;
}

void Stress::stres_gradcorr() 
{
     
     if (xcf.igcx == 0  &&  xcf.igcc == 0)
     {
        return;
     } 
     double sigma_gradcorr[3][3];
     double* p= &sigma_gradcorr[0][0];
     for(int i=0;i<9;i++)
        *p++ = 0;

     bool igcc_is_lyp = false;
     if( xcf.igcc == 3 || xcf.igcc == 7)
     {
                igcc_is_lyp = true;
     }

     assert(NSPIN>0);
     const double fac = 1.0/ NSPIN;

     // doing FFT to get rho in G space: rhog1 
     CHR.set_rhog(CHR.rho[0], CHR.rhog[0]);
     if(NSPIN==2)//mohan fix bug 2012-05-28
     {
          CHR.set_rhog(CHR.rho[1], CHR.rhog[1]);
     }
     CHR.set_rhog(CHR.rho_core, CHR.rhog_core);
    
	double* rhotmp1;
	double* rhotmp2;
	complex<double>* rhogsum1;
	complex<double>* rhogsum2;
	Vector3<double>* gdr1;
	Vector3<double>* gdr2;
 
     rhotmp1 = new double[pw.nrxx];
     rhogsum1 = new complex<double>[pw.ngmc];
     ZEROS(rhotmp1, pw.nrxx);
     ZEROS(rhogsum1, pw.ngmc);
     for(int ir=0; ir<pw.nrxx; ir++) rhotmp1[ir] = CHR.rho[0][ir] + fac * CHR.rho_core[ir];
     for(int ig=0; ig<pw.ngmc; ig++) rhogsum1[ig] = CHR.rhog[0][ig] + fac * CHR.rhog_core[ig];
	gdr1 = new Vector3<double>[pw.nrxx];
	ZEROS(gdr1, pw.nrxx);

        GGA_PW::grad_rho( rhogsum1 , gdr1 );

        if(NSPIN==2)
        {
                rhotmp2 = new double[pw.nrxx];
                rhogsum2 = new complex<double>[pw.ngmc];
                ZEROS(rhotmp2, pw.nrxx);
                ZEROS(rhogsum2, pw.ngmc);
                for(int ir=0; ir<pw.nrxx; ir++) rhotmp2[ir] = CHR.rho[1][ir] + fac * CHR.rho_core[ir];
                for(int ig=0; ig<pw.ngmc; ig++) rhogsum2[ig] = CHR.rhog[1][ig] + fac * CHR.rhog_core[ig];

                        gdr2 = new Vector3<double>[pw.nrxx];
                        ZEROS(gdr2, pw.nrxx);

                GGA_PW::grad_rho( rhogsum2 , gdr2 );
        }
        
        const double epsr = 1.0e-6;
        const double epsg = 1.0e-10;

        double grho2a = 0.0;
        double grho2b = 0.0;
        double sx = 0.0;
        double sc = 0.0;
        double v1x = 0.0;
        double v2x = 0.0;
        double v1c = 0.0;
        double v2c = 0.0;
        double vtxcgc = 0.0;
        double etxcgc = 0.0;

        if(NSPIN==1||NSPIN==4)
        {
                double segno;
                for(int ir=0; ir<pw.nrxx; ir++)
                {
                        const double arho = std::abs( rhotmp1[ir] );
			if(arho > epsr)
			{
				grho2a = gdr1[ir].norm2();
				if( grho2a > epsg )
				{
					if( rhotmp1[ir] >= 0.0 ) segno = 1.0;
					if( rhotmp1[ir] < 0.0 ) segno = -1.0;

					XC_Functional::gcxc( arho, grho2a, sx, sc, v1x, v2x, v1c, v2c);
					double tt[3];
					tt[0] = gdr1[ir].x;
					tt[1] = gdr1[ir].y;
					tt[2] = gdr1[ir].z;
					for(int l = 0;l< 3;l++){
					   for(int m = 0;m< l+1;m++){
						sigma_gradcorr[l][m] += tt[l] * tt[m] * e2 * (v2x + v2c);
					   }
					}
				}
			}
		} 
	}
	else if(NSPIN==2)
	{
		double v1cup = 0.0;
		double v1cdw = 0.0;
		double v2cup = 0.0;
		double v2cdw = 0.0;
		double v1xup = 0.0;
		double v1xdw = 0.0;
		double v2xup = 0.0;
		double v2xdw = 0.0;
		double v2cud = 0.0;
		double v2c = 0.0;
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			double rh = rhotmp1[ir] + rhotmp2[ir];
			grho2a = gdr1[ir].norm2();;
			grho2b = gdr2[ir].norm2();;
			//XC_Functional::gcx_spin();
			gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
				sx, v1xup, v1xdw, v2xup, v2xdw);

			if(rh > epsr)
			{
				if(igcc_is_lyp)
				{
					WARNING_QUIT("stress","igcc_is_lyp is not available now.");
				}
				else
				{
					double zeta = ( rhotmp1[ir] - rhotmp2[ir] ) / rh;
					double grh2 = (gdr1[ir]+gdr2[ir]).norm2();
					//XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					v2cup = v2c;
					v2cdw = v2c;
					v2cud = v2c;
				}
			}
			else
			{
				sc = 0.0;
				v1cup = 0.0;
				v1cdw = 0.0;
				v2c = 0.0;
				v2cup = 0.0;
				v2cdw = 0.0;
				v2cud = 0.0;
			}
			double tt1[3],tt2[3];
			{
				tt1[0] = gdr1[ir].x;
				tt1[1] = gdr1[ir].y;
				tt1[2] = gdr1[ir].z;
				tt2[0] = gdr2[ir].x;
				tt2[1] = gdr2[ir].y;
				tt2[2] = gdr2[ir].z;
			}
			for(int l = 0;l< 3;l++){
			    for(int m = 0;m< l+1;m++){
				//    exchange
				sigma_gradcorr [l][m] += tt1[l] * tt1[m] * e2 * v2xup + 
							tt2[l] * tt2[m] * e2 * v2xdw;
				//    correlation
				sigma_gradcorr [l][m] += ( tt1[l] * tt1[m] * v2cup + 
							tt2[l] * tt2[m] * v2cdw + 
							(tt1[l] * tt2[m] +
							tt2[l] * tt1[m] ) * v2cud ) * e2;
			    }
			}
		}
	}

	for(int l = 0;l< 3;l++){
		for(int m = 0;m< l;m++){
			sigma_gradcorr[m][l] = sigma_gradcorr[l][m];
		}
	}
	for(int l = 0;l<3;l++){
		for(int m = 0;m<3;m++){
			Parallel_Reduce::reduce_double_pool( sigma_gradcorr[l][m] );
		}
	}
	p= &sigma_gradcorr[0][0];
	double* p1 = &sigmaxc[0][0];
	for(int i=0;i<9;i++){
		*p /= pw.ncxyz ;
		*p1++ += *p++;  
	}

	delete[] rhotmp1;
	delete[] rhogsum1;
	delete[] gdr1;
	if(NSPIN==2)
	{
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
	}
	return;
}


void Stress::stres_cc()
{
    
    int nt,ng,l,m,ir;
    double fact;
    complex<double> sigmadiag;
    double* rhocg;
    double g[3];
    
       
    int judge=0;
    for(nt=0;nt<ucell.ntype;nt++){
       if(ucell.atoms[nt].nlcc) judge++;
    }
    if(judge==0) return;
           
    //recalculate the exchange-correlation potential
    matrix vxc(NSPIN, pw.nrxx);
    H_XC_pw::v_xc(pw.nrxx, pw.ncxyz, ucell.omega, CHR.rho, CHR.rho_core, vxc);

    complex<double> * psic = new complex<double> [pw.nrxx];
    ZEROS(psic, pw.nrxx);
 
    if(NSPIN==1||NSPIN==4){
       for(ir=0;ir<pw.nrxx;ir++){
         // psic[ir] = vxc(0,ir);
         psic[ir] = complex<double>(vxc(0, ir),  0.0);
       }
    }
    else {
       for(ir=0;ir<pw.nrxx;ir++){
          psic[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
       }
    }
          // to G space
    pw.FFT_chg.FFT3D(psic, -1);

    //psic cantains now Vxc(G)
    rhocg= new double [pw.nggm];
    ZEROS(rhocg, pw.nggm);

    sigmadiag=0.0;
    if(INPUT.gamma_only){
       fact=2.0;
    }
    else{
       fact=1.0;
    }
    for(nt=0;nt<ucell.ntype;nt++){
        if(ucell.atoms[nt].nlcc){
        //drhoc();
        CHR.non_linear_core_correction(
                ppcell.numeric,
                ucell.atoms[nt].msh,
                ucell.atoms[nt].r,
                ucell.atoms[nt].rab,
                ucell.atoms[nt].rho_atc,
                rhocg);
        
  
        //diagonal term 
        if (pw.gstart==1) sigmadiag += conj(psic [pw.ig2fftc[0]] ) * pw.strucFac (nt, 0) * rhocg [pw.ig2ngg[0] ];
        for( ng = pw.gstart;ng< pw.ngmc;ng++){
           sigmadiag +=  conj(psic[pw.ig2fftc[ng]] ) * 
                pw.strucFac (nt, ng) * rhocg [pw.ig2ngg[ng] ] * fact;
        }
        
        
        deriv_drhoc (
             ppcell.numeric,
             ucell.atoms[nt].msh, 
             ucell.atoms[nt].r, 
             ucell.atoms[nt].rab, 
             ucell.atoms[nt].rho_atc, 
             rhocg);
        // non diagonal term (g=0 contribution missing)
        for( ng = pw.gstart;ng< pw.ngmc;ng++){
           g[0] = pw.gcar[ng].x;
           g[1] = pw.gcar[ng].y;
           g[2] = pw.gcar[ng].z;
           for( l = 0;l< 3;l++){
              for (m = 0;m< 3;m++){
                 const complex<double> t = conj(psic[pw.ig2fftc[ng]] )
                      * pw.strucFac (nt, ng) * rhocg [pw.ig2ngg[ng] ] * ucell.tpiba *
                      g [l] * g [m] / pw.gcar[ng].norm() * fact;
                 sigmaxcc [l][ m] += t.real();
              }//end m
           }//end l
        }//end ng
     }//end if
   }//end nt
   for( l = 0;l< 3;l++){
     sigmaxcc [l][ l] += sigmadiag.real();
   }
   //mp_sum(  sigmaxcc, intra_bgrp_comm );
   for( l = 0;l< 3;l++){
              for (m = 0;m< 3;m++){
                    Parallel_Reduce::reduce_double_pool( sigmaxcc[l][m] );
              }
   }

   delete[] rhocg;
   delete[] psic;
   return;
}


void Stress::deriv_drhoc (
    const bool &numeric,
    const int mesh,
    const double *r,
    const double *rab,
    const double *rhoc,
    double *drhocg)
{
  //int ngl, mesh;
  // input: the number of g shell
  // input: the number of radial mesh points

  // input: the number of G shells
  // input: the radial mesh
  // input: the derivative of the radial mesh
  // input: the radial core charge
  // output: fourier transform of d Rho_c/dG
  //
  //     here the local variables
  //
  double gx = 0, rhocg1 = 0;
  // the modulus of g for a given shell
  // the fourier transform
  double *aux;
  // auxiliary memory for integration

  int  ir, igl, igl0;
  // counter on radial mesh points
  // counter on g shells
  // lower limit for loop on ngl

  //
  // G=0 term
  //
  if (pw.ggs[0] < 1.0e-8){
     drhocg [0] = 0.0;
     igl0 = 1;
  }
  else{
     igl0 = 0;
  }
  //
  // G <> 0 term
  //
  aux= new double[ mesh];

  for( igl = igl0;igl< pw.nggm;igl++){
     gx = sqrt(pw.ggs [igl] * ucell.tpiba2);
     for( ir = 0;ir< mesh; ir++){
        aux [ir] = r [ir] * rhoc [ir] * (r [ir] * cos (gx * r [ir] ) / gx - sin (gx * r [ir] ) / pow(gx,2));
     }//ir
     Mathzone::Simpson_Integral(mesh, aux, rab, rhocg1);
     drhocg [igl] = FOUR_PI / ucell.omega * rhocg1;
  }//igl
  delete [] aux;

  return;
}

