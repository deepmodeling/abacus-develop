#include "./stress_func.h"
#include "./myfunc.h"
#include "./H_Hartree_pw.h"

//calculate the Hartree part in PW or LCAO base
void Stress_Func::stress_har(matrix& sigma, const bool is_pw)
{
	timer::tick("Stress_Func","stress_har",'F');
	double shart,g2;
	const double eps=1e-8;
	int l,m,nspin0;

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
		daxpy (pw.nrxx, 1.0, CHR.rho[is], 1, psic0, 2);
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
//	double ehart=0;
	for (int ig = pw.gstart; ig<pw.ngmc; ig++)
	{
		const int j = pw.ig2fftc[ig];
		const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);

//		ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
	//            vh_g[ig] = fac * Porter[j];
		shart= ( conj( Porter[j] ) * Porter[j] ).real()/(ucell.tpiba2 * pw.gg [ig]);
        //test

        //   cout<<"g "<<g[0]<<" "<<g[1]<<" "<<g[2]<<endl;
        //   cout<<"gg "<<pw.gg[i]<<" "<<pw.gcar[i].norm2()<<endl;
        //   cout<<"Porter "<<Porter[j]<<endl;
        //   cout<<"tpiba2 "<<ucell.tpiba2<<endl;


		for(int l=0;l<3;l++)
		{
			for(int m=0;m<l+1;m++)
			{
				sigma(l, m) += shart * 2 * pw.get_G_cartesian_projection(ig, l) * pw.get_G_cartesian_projection(ig, m) / pw.gg[ig];
			}
		}
	}
//	Parallel_Reduce::reduce_double_pool( en.ehart );
//	ehart *= 0.5 * ucell.omega;
        //cout<<"ehart "<<ehart<<" en.ehart "<< en.ehart<<endl;
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<l+1;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}

//        Parallel_Reduce::reduce_double_pool( ehart );
//        ehart *= 0.5 * ucell.omega;
        //psic(:)=(0.0,0.0)
	if(is_pw&&INPUT.gamma_only)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				sigma(l,m) *= e2 * FOUR_PI;
			}
		}
	}
	else
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				sigma(l,m) *= 0.5 * e2 * FOUR_PI;
			}
		}
	}
	
	for(int l=0;l<3;l++)
	{
		if(is_pw) sigma(l,l) -= H_Hartree_pw::hartree_energy /ucell.omega;
		else sigma(l,l) += H_Hartree_pw::hartree_energy /ucell.omega;
		for(int m=0;m<l;m++)
		{
			sigma(m,l)=sigma(l,m);
		}
	}
	
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigma(l,m)*=-1;
		}
	}

	delete[] vh_g;
	delete[] psic;
	delete[] psic0;
	timer::tick("Stress_Func","stress_har");
	return;
}
