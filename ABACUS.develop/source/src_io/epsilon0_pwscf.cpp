//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-30
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating the macroscopic dielectric tensor  using the linear
//response theory .
//-------------------------------------------------------------------------------
#include "src_pw/global.h"
#include "src_pw/hamilt_pw.h"
#include "../src_pw/wavefunc_in_pw.h"
#include "src_io/optical.h"
#include "epsilon0_pwscf.h"
#include <iostream>
#include <math.h>

using namespace std;

Epsilon0_pwscf epsilon0_pwscf;

Epsilon0_pwscf::Epsilon0_pwscf()
          :init_finish(false)
{
      epsilon = false;
	  intersmear = 0.01;	  
	  intrasmear = 0.0;
	  domega = 0.01;
	  nomega = 300;
	  shift = 0.0;
	  metalcalc = false;
	  degauss = 0.01;			  
}
		  
Epsilon0_pwscf::~Epsilon0_pwscf()
{	
}

void Epsilon0_pwscf:: Cal_epsilon0()
{
	cout << "intersmear = " << intersmear << endl;
	cout << "intrasmear = " << intrasmear << endl;
	cout << "domega = "<<domega<<endl;
	cout << "nomega = "<<nomega<<endl;
	cout << "shift = "<< shift << endl;
	cout << "metalcalc = "<<metalcalc<<endl;
	cout << "degauss = "<<degauss<<endl;
	
	if( !init_finish )
	{
		Init();
		init_finish = true;
	}
	
	for(int i=0; i<9; i++)
		for(int j=0; j<nomega; j++)
		{
			epsr[i][j] = 0.0;
			epsi[i][j] = 0.0;
		}

		
	if(NSPIN == 1)
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			Cal_dipole(ik);
			for(int ib1=0; ib1<NBANDS; ib1++)
				for(int ib2=0; ib2<NBANDS; ib2++)
				{
					dipole[0][ib1][ib2] = dipole_aux[0][ib1][ib2] * conj(dipole_aux[0][ib1][ib2]);
					dipole[1][ib1][ib2] = dipole_aux[0][ib1][ib2] * conj(dipole_aux[1][ib1][ib2]);
					dipole[2][ib1][ib2] = dipole_aux[0][ib1][ib2] * conj(dipole_aux[2][ib1][ib2]);
					dipole[3][ib1][ib2] = dipole_aux[1][ib1][ib2] * conj(dipole_aux[0][ib1][ib2]);
					dipole[4][ib1][ib2] = dipole_aux[1][ib1][ib2] * conj(dipole_aux[1][ib1][ib2]);
					dipole[5][ib1][ib2] = dipole_aux[1][ib1][ib2] * conj(dipole_aux[2][ib1][ib2]);
					dipole[6][ib1][ib2] = dipole_aux[2][ib1][ib2] * conj(dipole_aux[0][ib1][ib2]);
					dipole[7][ib1][ib2] = dipole_aux[2][ib1][ib2] * conj(dipole_aux[1][ib1][ib2]);
					dipole[8][ib1][ib2] = dipole_aux[2][ib1][ib2] * conj(dipole_aux[2][ib1][ib2]);
				}
					
			for(int ib2=0; ib2<NBANDS; ib2++)
			{
				//cout <<"ik= "<<ik<<" ib2= "<<ib2<<" focc= "<<focc(ib2,ik)<<endl;
				if(focc(ib2,ik) < 2.0)
				{
					for(int ib1=0; ib1<NBANDS; ib1++)
					{
						if(ib1 == ib2)
							continue;
						if(focc(ib1,ik) > 0.0001)
						{
							if( fabs(focc(ib2,ik)-focc(ib1,ik)) > 0.001)
							{
								double etrans = ( wf.ekb[ik][ib2] - wf.ekb[ik][ib1] ) + shift;
								for(int iw=0; iw<nomega; iw++)
								{
									double w = wgrid(iw);
									for(int j=0; j<9; j++)
									{
										epsi[j][iw] += dipole[j][ib1][ib2].real() * intersmear * w * focc(ib1,ik) /((pow(etrans * etrans - w * w,2) + intersmear * intersmear * w * w) * etrans);
										epsr[j][iw] += dipole[j][ib1][ib2].real() * focc(ib1,ik) * (etrans * etrans - w * w)/((pow(etrans * etrans - w * w,2) + intersmear * intersmear * w * w) * etrans);
									}

								}
							}
						}
					}
				}
			}
			
			if(metalcalc)
			{
				for(int ib1=0; ib1<NBANDS; ib1++)
				{
					if(focc(ib1,ik) < 2.0)
					{
						if(focc(ib1,ik) > 0.0001)
						{
							for(int iw=0; iw<nomega; iw++)
							{
								double w = wgrid(iw);
								for(int j=0; j<9; j++)
								{
									epsi[j][iw] += dipole[j][ib1][ib1].real() * intrasmear * w * exp((wf.ekb[ik][ib1]- en.ef)/degauss)/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/degauss/(w*w*w*w + intrasmear * intrasmear * w * w);
									epsr[j][iw] -= dipole[j][ib1][ib1].real() * w * w  * exp((wf.ekb[ik][ib1]- en.ef)/degauss)/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/degauss/(w*w*w*w + intrasmear * intrasmear * w * w);
								}
							}
						}
					}
				}
			}
		}
		
		double coff = 64 * PI /ucell.omega/kv.nks;
		cout << "all finish" << endl;
		
		for(int iw=0; iw<nomega; iw++)
			for(int j=0; j<9; j++)
			{
				epsi[j][iw] = epsi[j][iw] * coff;
				epsr[j][iw] = epsr[j][iw] * coff; 
			}
		
		for(int iw=0; iw<nomega; iw++)
			for(int j=0; j<9; j=j+4)
			{
				epsr[j][iw] += 1.0;
			}
			
		ofs_running<<endl;
		ofs_running<<" The real part of the macroscopic dielectric constant:"<<endl;
		ofs_running<<setw(15)<<"omega"<<setw(15)<<"XX"<<setw(15)<<"XY"<<setw(15)<<"XZ"<<setw(15)<<"YX"<<setw(15)<<"YY"<<setw(15)<<"YZ"<<setw(15)<<"ZX"<<setw(15)<<"ZY"<<setw(15)<<"ZZ"<<endl; 
		for(int i=0; i<nomega; i++)
		{
			ofs_running<<setprecision(2)<<setw(15)<<(i*domega)<<setprecision(2)<<setw(15)<<epsr[0][i]<<setprecision(2)<<setw(15)<<epsr[1][i]<<setprecision(2)<<setw(15)<<epsr[2][i]<<setprecision(2)<<setw(15)<<epsr[3][i]<<setprecision(2)<<setw(15)<<epsr[4][i]<<setprecision(2)<<setw(15)<<epsr[5][i]<<setprecision(2)<<setw(15)<<epsr[6][i]<<setprecision(2)<<setw(15)<<epsr[7][i]<<setprecision(2)<<setw(15)<<epsr[8][i]<<endl;
		}

		ofs_running<<endl; 
		ofs_running<<" The imag part of the macroscopic dielectric constant:"<<endl;
		ofs_running<<setw(15)<<"omega"<<setw(15)<<"XX"<<setw(15)<<"XY"<<setw(15)<<"XZ"<<setw(15)<<"YX"<<setw(15)<<"YY"<<setw(15)<<"YZ"<<setw(15)<<"ZX"<<setw(15)<<"ZY"<<setw(15)<<"ZZ"<<endl;  
		for(int i=0; i<nomega; i++)
		{
			ofs_running<<setprecision(2)<<setw(15)<<(i*domega)<<setprecision(2)<<setw(15)<<epsi[0][i]<<setprecision(2)<<setw(15)<<epsi[1][i]<<setprecision(2)<<setw(15)<<epsi[2][i]<<setprecision(2)<<setw(15)<<epsi[3][i]<<setprecision(2)<<setw(15)<<epsi[4][i]<<setprecision(2)<<setw(15)<<epsi[5][i]<<setprecision(2)<<setw(15)<<epsi[6][i]<<setprecision(2)<<setw(15)<<epsi[7][i]<<setprecision(2)<<setw(15)<<epsi[8][i]<<endl;
		}					
		
		/*for(int iw=0; iw<nomega; iw++)
		{
			cout <<"epsi[0]["<<iw<<"] = "<< epsi[0][iw]<<endl;
			cout <<"epsr[0]["<<iw<<"] = "<< epsr[0][iw]<<endl;
		}*/

		ofs_running<<" Macroscopic dielectric constant matrix :"<<endl;
		ofs_running<<setprecision(2)<<setw(15)<<epsr[0][0]<<setprecision(2)<<setw(15)<<epsr[1][0]<<setprecision(2)<<setw(15)<<epsr[2][0]<<endl;
		ofs_running<<setprecision(2)<<setw(15)<<epsr[3][0]<<setprecision(2)<<setw(15)<<epsr[4][0]<<setprecision(2)<<setw(15)<<epsr[5][0]<<endl;
		ofs_running<<setprecision(2)<<setw(15)<<epsr[6][0]<<setprecision(2)<<setw(15)<<epsr[7][0]<<setprecision(2)<<setw(15)<<epsr[8][0]<<endl;	
		
		//cout <<"Macroscopic dielectric constant matrix :"<<endl;
		//cout << epsr[0][0] <<"  "<<epsr[1][0] <<"  "<<epsr[2][0]<<"  "<<endl;
		//cout << epsr[3][0] <<"  "<<epsr[4][0] <<"  "<<epsr[5][0]<<"  "<<endl;
		//cout << epsr[6][0] <<"  "<<epsr[7][0] <<"  "<<epsr[8][0]<<"  "<<endl;
			
	}
	else if(NSPIN == 2)
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			Cal_dipole(ik);
			for(int ib1=0; ib1<NBANDS; ib1++)
				for(int ib2=0; ib2<NBANDS; ib2++)
				{
					dipole[0][ib1][ib2] = dipole_aux[0][ib1][ib2] * conj(dipole_aux[0][ib1][ib2]);
					dipole[1][ib1][ib2] = dipole_aux[0][ib1][ib2] * conj(dipole_aux[1][ib1][ib2]);
					dipole[2][ib1][ib2] = dipole_aux[0][ib1][ib2] * conj(dipole_aux[2][ib1][ib2]);
					dipole[3][ib1][ib2] = dipole_aux[1][ib1][ib2] * conj(dipole_aux[0][ib1][ib2]);
					dipole[4][ib1][ib2] = dipole_aux[1][ib1][ib2] * conj(dipole_aux[1][ib1][ib2]);
					dipole[5][ib1][ib2] = dipole_aux[1][ib1][ib2] * conj(dipole_aux[2][ib1][ib2]);
					dipole[6][ib1][ib2] = dipole_aux[2][ib1][ib2] * conj(dipole_aux[0][ib1][ib2]);
					dipole[7][ib1][ib2] = dipole_aux[2][ib1][ib2] * conj(dipole_aux[1][ib1][ib2]);
					dipole[8][ib1][ib2] = dipole_aux[2][ib1][ib2] * conj(dipole_aux[2][ib1][ib2]);
				}
					
			for(int ib2=0; ib2<NBANDS; ib2++)
			{
				//cout <<"ik= "<<ik<<" ib2= "<<ib2<<" focc= "<<focc(ib2,ik)<<endl;
				if(focc(ib2,ik) < 1.0)
				{
					for(int ib1=0; ib1<NBANDS; ib1++)
					{
						if(ib1 == ib2)
							continue;
						if(focc(ib1,ik) > 0.0001)
						{
							if( fabs(focc(ib2,ik)-focc(ib1,ik)) > 0.001)
							{
								double etrans = ( wf.ekb[ik][ib2] - wf.ekb[ik][ib1] ) + shift;
								for(int iw=0; iw<nomega; iw++)
								{
									double w = wgrid(iw);
									for(int j=0; j<9; j++)
									{
										epsi[j][iw] += dipole[j][ib1][ib2].real() * intersmear * w * focc(ib1,ik) /((pow(etrans * etrans - w * w,2) + intersmear * intersmear * w * w) * etrans);
										epsr[j][iw] += dipole[j][ib1][ib2].real() * focc(ib1,ik) * (etrans * etrans - w * w)/((pow(etrans * etrans - w * w,2) + intersmear * intersmear * w * w) * etrans);
									}

								}
							}
						}
					}
				}
			}
			
			if(metalcalc)
			{
				for(int ib1=0; ib1<NBANDS; ib1++)
				{
					if(focc(ib1,ik) < 1.0)
					{
						if(focc(ib1,ik) > 0.0001)
						{
							for(int iw=0; iw<nomega; iw++)
							{
								double w = wgrid(iw);
								for(int j=0; j<9; j++)
								{
									epsi[j][iw] += dipole[j][ib1][ib1].real() * intrasmear * w * exp((wf.ekb[ik][ib1]- en.ef)/degauss)/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/degauss/(w*w*w*w + intrasmear * intrasmear * w * w);
									epsr[j][iw] -= dipole[j][ib1][ib1].real() * w * w  * exp((wf.ekb[ik][ib1]- en.ef)/degauss)/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/(1.0 + exp((wf.ekb[ik][ib1]- en.ef)/degauss))/degauss/(w*w*w*w + intrasmear * intrasmear * w * w);
								}
							}
						}
					}
				}
			}
		}
		
		double coff = 128 * PI /ucell.omega/kv.nks;
		cout << "all finish" << endl;

		for(int iw=0; iw<nomega; iw++)
			for(int j=0; j<9; j++)
			{
				epsi[j][iw] = epsi[j][iw] * coff;
				epsr[j][iw] = epsr[j][iw] * coff; 
			}
		
		for(int iw=0; iw<nomega; iw++)
			for(int j=0; j<9; j=j+4)
			{
				epsr[j][iw] += 1.0;
			}
			
		ofs_running<<endl;
		ofs_running<<" The real part of the macroscopic dielectric constant:"<<endl;
		ofs_running<<setw(15)<<"omega"<<setw(15)<<"XX"<<setw(15)<<"XY"<<setw(15)<<"XZ"<<setw(15)<<"YX"<<setw(15)<<"YY"<<setw(15)<<"YZ"<<setw(15)<<"ZX"<<setw(15)<<"ZY"<<setw(15)<<"ZZ"<<endl; 
		for(int i=0; i<nomega; i++)
		{
			ofs_running<<setprecision(2)<<setw(15)<<(i*domega)<<setprecision(2)<<setw(15)<<epsr[0][i]<<setprecision(2)<<setw(15)<<epsr[1][i]<<setprecision(2)<<setw(15)<<epsr[2][i]<<setprecision(2)<<setw(15)<<epsr[3][i]<<setprecision(2)<<setw(15)<<epsr[4][i]<<setprecision(2)<<setw(15)<<epsr[5][i]<<setprecision(2)<<setw(15)<<epsr[6][i]<<setprecision(2)<<setw(15)<<epsr[7][i]<<setprecision(2)<<setw(15)<<epsr[8][i]<<endl;
		}		

		ofs_running<<endl; 
		ofs_running<<" The imag part of the macroscopic dielectric constant:"<<endl;
		ofs_running<<setw(15)<<"omega"<<setw(15)<<"XX"<<setw(15)<<"XY"<<setw(15)<<"XZ"<<setw(15)<<"YX"<<setw(15)<<"YY"<<setw(15)<<"YZ"<<setw(15)<<"ZX"<<setw(15)<<"ZY"<<setw(15)<<"ZZ"<<endl;  
		for(int i=0; i<nomega; i++)
		{
			ofs_running<<setprecision(2)<<setw(15)<<(i*domega)<<setprecision(2)<<setw(15)<<epsi[0][i]<<setprecision(2)<<setw(15)<<epsi[1][i]<<setprecision(2)<<setw(15)<<epsi[2][i]<<setprecision(2)<<setw(15)<<epsi[3][i]<<setprecision(2)<<setw(15)<<epsi[4][i]<<setprecision(2)<<setw(15)<<epsi[5][i]<<setprecision(2)<<setw(15)<<epsi[6][i]<<setprecision(2)<<setw(15)<<epsi[7][i]<<setprecision(2)<<setw(15)<<epsi[8][i]<<endl;
		}		
		
		/*for(int iw=0; iw<nomega; iw++)
		{
			cout <<"epsi[0]["<<iw<<"] = "<< epsi[0][iw]<<endl;
			cout <<"epsr[0]["<<iw<<"] = "<< epsr[0][iw]<<endl;
		}*/

		ofs_running<<" Macroscopic dielectric constant matrix :"<<endl;
		ofs_running<<setprecision(2)<<setw(15)<<epsr[0][0]<<setprecision(2)<<setw(15)<<epsr[1][0]<<setprecision(2)<<setw(15)<<epsr[2][0]<<endl;
		ofs_running<<setprecision(2)<<setw(15)<<epsr[3][0]<<setprecision(2)<<setw(15)<<epsr[4][0]<<setprecision(2)<<setw(15)<<epsr[5][0]<<endl;
		ofs_running<<setprecision(2)<<setw(15)<<epsr[6][0]<<setprecision(2)<<setw(15)<<epsr[7][0]<<setprecision(2)<<setw(15)<<epsr[8][0]<<endl;	
			
		//cout <<"Macroscopic dielectric constant matrix :"<<endl;
		//cout << epsr[0][0] <<"  "<<epsr[1][0] <<"  "<<epsr[2][0]<<"  "<<endl;
		//cout << epsr[3][0] <<"  "<<epsr[4][0] <<"  "<<epsr[5][0]<<"  "<<endl;
		//cout << epsr[6][0] <<"  "<<epsr[7][0] <<"  "<<epsr[8][0]<<"  "<<endl;		
	}
	
	Delete();
	
	return;	
}

void Epsilon0_pwscf:: Init()
{
    
	dipole_aux = new complex<double>**[3];
	for(int i=0; i<3; i++)
	{
		dipole_aux[i] = new complex<double>*[NBANDS];
		for(int j=0; j<NBANDS; j++)
		{
			dipole_aux[i][j] = new complex<double>[NBANDS];
		}
	}
	
	dipole = new complex<double>**[9];
	for(int i=0; i<9; i++)
	{
		dipole[i] = new complex<double>*[NBANDS];
		for(int j=0; j<NBANDS; j++)
		{
			dipole[i][j] = new complex<double>[NBANDS];
		}
	}
	
	epsi = new double*[9];
	for(int i=0; i<9; i++)
	{
		epsi[i] = new double[nomega];
	}
	
	epsr = new double*[9];
	for(int i=0; i<9; i++)
	{
		epsr[i] = new double[nomega];
	}
	
	return;
}

void Epsilon0_pwscf:: Delete()
{
	
	if(init_finish)
	{
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<NBANDS; j++)
			{
				delete[] dipole_aux[i][j];
			}
			delete[] dipole_aux[i];
		}
		delete[] dipole_aux;
		
		for(int i=0; i<9; i++)
		{
			for(int j=0; j<NBANDS; j++)
			{
				delete[] dipole[i][j];
			}
			delete[] dipole[i];
		}
		delete[] dipole;
		
		for(int i=0; i<9; i++)
		{
			delete[] epsr[i];
		}
		delete[] epsr;
		
		for(int i=0; i<9; i++)
		{
			delete[] epsi[i];
		}
		delete[] epsi;
	}
}


double Epsilon0_pwscf:: wgrid(int iw)
{
	return (iw * domega);
}

double Epsilon0_pwscf:: focc(int ib, int ik)
{
	if(NSPIN == 1)
		return (wf.wg(ik,ib) * 2.0/ kv.wk[ik] );
	else if(NSPIN == 2)
		return (wf.wg(ik,ib) * 1.0/ kv.wk[ik] );
	else				// Peize Lin add 2019-05-01
		throw domain_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
} 
void Epsilon0_pwscf:: Cal_dipole(int ik)
{
	
	complex<double> dipole_aux_core[3][NBANDS][NBANDS];
	
	for(int i=0; i<3; i++)
		for(int ib1=0; ib1<NBANDS; ib1++)
			for(int ib2=0; ib2<NBANDS; ib2++)
			{
				dipole_aux_core[i][ib1][ib2] = complex<double>(0.0,0.0);
			}
	
	
	for(int ib2=0; ib2<NBANDS; ib2++)
	{
		if( focc(ib2,ik) < 2.0)
		{
			for(int ib1=0; ib1<NBANDS; ib1++)
			{
				if(ib1 == ib2)
					continue;
				if( focc(ib1,ik) > 0.0001)
				{
					for(int ig=0; ig<kv.ngk[ik]; ig++)
					{
						dipole_aux_core[0][ib1][ib2] += conj(wf.evc[ik](ib1, ig)) * (pw.get_G_cartesian_projection(ig, 0)) * (TWO_PI / ucell.lat0)) * wf.evc[ik](ib2, ig);
						dipole_aux_core[1][ib1][ib2] += conj(wf.evc[ik](ib1,ig)) * (pw.get_G_cartesian_projection(ig, 1))*(TWO_PI/ucell.lat0)) * wf.evc[ik](ib2,ig);
						dipole_aux_core[2][ib1][ib2] += conj(wf.evc[ik](ib1,ig)) * (pw.get_G_cartesian_projection(ig, 2))*(TWO_PI/ucell.lat0)) * wf.evc[ik](ib2,ig);
					}
				}
			}
			
		}
	}
	
	if(metalcalc)
	{
		for(int ib1=0; ib1<NBANDS; ib1++)
		{
			for(int ig=0; ig<kv.ngk[ik]; ig++)
			{
				dipole_aux_core[0][ib1][ib1] += conj(wf.evc[ik](ib1, ig)) * (pw.get_GPlusK_cartesian_projection(ik, ig, 0) * (TWO_PI / ucell.lat0) * wf.evc[ik](ib1,ig);
				dipole_aux_core[1][ib1][ib1] += conj(wf.evc[ik](ib1, ig)) * (pw.get_GPlusK_cartesian_projection(ik, ig, 1) * (TWO_PI / ucell.lat0) * wf.evc[ik](ib1, ig);
				dipole_aux_core[2][ib1][ib1] += conj(wf.evc[ik](ib1, ig)) * (pw.get_GPlusK_cartesian_projection(ik, ig, 2) * (TWO_PI / ucell.lat0) * wf.evc[ik](ib1, ig);
			}
		}
	}
	
	double dipole_aux_core_R[3][NBANDS][NBANDS];
	double dipole_aux_core_I[3][NBANDS][NBANDS];
	double dipole_aux_R[3][NBANDS][NBANDS];
	double dipole_aux_I[3][NBANDS][NBANDS];
	
	for(int i=0; i<3; i++)
		for(int ib1=0; ib1<NBANDS; ib1++)
			for(int ib2=0; ib2<NBANDS; ib2++)
			{
				dipole_aux_core_R[i][ib1][ib2] = dipole_aux_core[i][ib1][ib2].real();
				dipole_aux_core_I[i][ib1][ib2] = dipole_aux_core[i][ib1][ib2].imag();
			}
			
#ifdef __MPI
	MPI_Allreduce(dipole_aux_core_R,dipole_aux_R,3*NBANDS*NBANDS,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	MPI_Allreduce(dipole_aux_core_I,dipole_aux_I,3*NBANDS*NBANDS,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif
	
	for(int i=0; i<3; i++)
		for(int ib1=0; ib1<NBANDS; ib1++)
			for(int ib2=0; ib2<NBANDS; ib2++)
			{
				dipole_aux[i][ib1][ib2] = complex<double>( dipole_aux_R[i][ib1][ib2], dipole_aux_I[i][ib1][ib2]);
			}
			
	
	return;
}
