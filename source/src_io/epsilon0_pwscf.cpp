//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-30
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating the macroscopic dielectric tensor  using the linear
//response theory .
//-------------------------------------------------------------------------------
#include "../src_pw/global.h"
#include "../src_pw/hamilt.h"
#include "../src_lcao/wavefunc_in_pw.h"
#include "optical.h"
#include "epsilon0_pwscf.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace GlobalC
{
Epsilon0_pwscf epsilon0_pwscf;
}

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
	  smearing_sigma = 0.01;
}

Epsilon0_pwscf::~Epsilon0_pwscf()
{
}

void Epsilon0_pwscf:: Cal_epsilon0()
{
	std::cout << "intersmear = " << intersmear << std::endl;
	std::cout << "intrasmear = " << intrasmear << std::endl;
	std::cout << "domega = "<<domega<<std::endl;
	std::cout << "nomega = "<<nomega<<std::endl;
	std::cout << "shift = "<< shift << std::endl;
	std::cout << "metalcalc = "<<metalcalc<<std::endl;
	std::cout << "smearing_sigma = "<<smearing_sigma<<std::endl;

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


	if(GlobalV::NSPIN == 1)
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			Cal_dipole(ik);
			for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
				for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
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

			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				//std::cout <<"ik= "<<ik<<" ib2= "<<ib2<<" focc= "<<focc(ib2,ik)<<std::endl;
				if(focc(ib2,ik) < 2.0)
				{
					for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
					{
						if(ib1 == ib2)
							continue;
						if(focc(ib1,ik) > 0.0001)
						{
							if( fabs(focc(ib2,ik)-focc(ib1,ik)) > 0.001)
							{
								double etrans = ( GlobalC::wf.ekb[ik][ib2] - GlobalC::wf.ekb[ik][ib1] ) + shift;
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
				for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
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
									epsi[j][iw] += dipole[j][ib1][ib1].real() * intrasmear * w * exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma)/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/smearing_sigma/(w*w*w*w + intrasmear * intrasmear * w * w);
									epsr[j][iw] -= dipole[j][ib1][ib1].real() * w * w  * exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma)/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/smearing_sigma/(w*w*w*w + intrasmear * intrasmear * w * w);
								}
							}
						}
					}
				}
			}
		}
		
		double coff = 64 * ModuleBase::PI /GlobalC::ucell.omega/GlobalC::kv.nks;
		std::cout << "all finish" << std::endl;

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

		GlobalV::ofs_running<<std::endl;
		GlobalV::ofs_running<<" The real part of the macroscopic dielectric constant:"<<std::endl;
		GlobalV::ofs_running<<std::setw(15)<<"omega"<<std::setw(15)<<"XX"<<std::setw(15)<<"XY"<<std::setw(15)<<"XZ"<<std::setw(15)<<"YX"<<std::setw(15)<<"YY"<<std::setw(15)<<"YZ"<<std::setw(15)<<"ZX"<<std::setw(15)<<"ZY"<<std::setw(15)<<"ZZ"<<std::endl;
		for(int i=0; i<nomega; i++)
		{
			GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<(i*domega)<<std::setprecision(2)<<std::setw(15)<<epsr[0][i]<<std::setprecision(2)<<std::setw(15)<<epsr[1][i]<<std::setprecision(2)<<std::setw(15)<<epsr[2][i]<<std::setprecision(2)<<std::setw(15)<<epsr[3][i]<<std::setprecision(2)<<std::setw(15)<<epsr[4][i]<<std::setprecision(2)<<std::setw(15)<<epsr[5][i]<<std::setprecision(2)<<std::setw(15)<<epsr[6][i]<<std::setprecision(2)<<std::setw(15)<<epsr[7][i]<<std::setprecision(2)<<std::setw(15)<<epsr[8][i]<<std::endl;
		}

		GlobalV::ofs_running<<std::endl;
		GlobalV::ofs_running<<" The imag part of the macroscopic dielectric constant:"<<std::endl;
		GlobalV::ofs_running<<std::setw(15)<<"omega"<<std::setw(15)<<"XX"<<std::setw(15)<<"XY"<<std::setw(15)<<"XZ"<<std::setw(15)<<"YX"<<std::setw(15)<<"YY"<<std::setw(15)<<"YZ"<<std::setw(15)<<"ZX"<<std::setw(15)<<"ZY"<<std::setw(15)<<"ZZ"<<std::endl;
		for(int i=0; i<nomega; i++)
		{
			GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<(i*domega)<<std::setprecision(2)<<std::setw(15)<<epsi[0][i]<<std::setprecision(2)<<std::setw(15)<<epsi[1][i]<<std::setprecision(2)<<std::setw(15)<<epsi[2][i]<<std::setprecision(2)<<std::setw(15)<<epsi[3][i]<<std::setprecision(2)<<std::setw(15)<<epsi[4][i]<<std::setprecision(2)<<std::setw(15)<<epsi[5][i]<<std::setprecision(2)<<std::setw(15)<<epsi[6][i]<<std::setprecision(2)<<std::setw(15)<<epsi[7][i]<<std::setprecision(2)<<std::setw(15)<<epsi[8][i]<<std::endl;
		}

		/*for(int iw=0; iw<nomega; iw++)
		{
			std::cout <<"epsi[0]["<<iw<<"] = "<< epsi[0][iw]<<std::endl;
			std::cout <<"epsr[0]["<<iw<<"] = "<< epsr[0][iw]<<std::endl;
		}*/

		GlobalV::ofs_running<<" Macroscopic dielectric constant matrix :"<<std::endl;
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<epsr[0][0]<<std::setprecision(2)<<std::setw(15)<<epsr[1][0]<<std::setprecision(2)<<std::setw(15)<<epsr[2][0]<<std::endl;
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<epsr[3][0]<<std::setprecision(2)<<std::setw(15)<<epsr[4][0]<<std::setprecision(2)<<std::setw(15)<<epsr[5][0]<<std::endl;
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<epsr[6][0]<<std::setprecision(2)<<std::setw(15)<<epsr[7][0]<<std::setprecision(2)<<std::setw(15)<<epsr[8][0]<<std::endl;

		//std::cout <<"Macroscopic dielectric constant matrix :"<<std::endl;
		//std::cout << epsr[0][0] <<"  "<<epsr[1][0] <<"  "<<epsr[2][0]<<"  "<<std::endl;
		//std::cout << epsr[3][0] <<"  "<<epsr[4][0] <<"  "<<epsr[5][0]<<"  "<<std::endl;
		//std::cout << epsr[6][0] <<"  "<<epsr[7][0] <<"  "<<epsr[8][0]<<"  "<<std::endl;

	}
	else if(GlobalV::NSPIN == 2)
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			Cal_dipole(ik);
			for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
				for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
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

			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				//std::cout <<"ik= "<<ik<<" ib2= "<<ib2<<" focc= "<<focc(ib2,ik)<<std::endl;
				if(focc(ib2,ik) < 1.0)
				{
					for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
					{
						if(ib1 == ib2)
							continue;
						if(focc(ib1,ik) > 0.0001)
						{
							if( fabs(focc(ib2,ik)-focc(ib1,ik)) > 0.001)
							{
								double etrans = ( GlobalC::wf.ekb[ik][ib2] - GlobalC::wf.ekb[ik][ib1] ) + shift;
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
				for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
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
									epsi[j][iw] += dipole[j][ib1][ib1].real() * intrasmear * w * exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma)/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/smearing_sigma/(w*w*w*w + intrasmear * intrasmear * w * w);
									epsr[j][iw] -= dipole[j][ib1][ib1].real() * w * w  * exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma)/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/(1.0 + exp((GlobalC::wf.ekb[ik][ib1]- GlobalC::en.ef)/smearing_sigma))/smearing_sigma/(w*w*w*w + intrasmear * intrasmear * w * w);
								}
							}
						}
					}
				}
			}
		}
		
		double coff = 128 * ModuleBase::PI /GlobalC::ucell.omega/GlobalC::kv.nks;
		std::cout << "all finish" << std::endl;

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

		GlobalV::ofs_running<<std::endl;
		GlobalV::ofs_running<<" The real part of the macroscopic dielectric constant:"<<std::endl;
		GlobalV::ofs_running<<std::setw(15)<<"omega"<<std::setw(15)<<"XX"<<std::setw(15)<<"XY"<<std::setw(15)<<"XZ"<<std::setw(15)<<"YX"<<std::setw(15)<<"YY"<<std::setw(15)<<"YZ"<<std::setw(15)<<"ZX"<<std::setw(15)<<"ZY"<<std::setw(15)<<"ZZ"<<std::endl;
		for(int i=0; i<nomega; i++)
		{
			GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<(i*domega)<<std::setprecision(2)<<std::setw(15)<<epsr[0][i]<<std::setprecision(2)<<std::setw(15)<<epsr[1][i]<<std::setprecision(2)<<std::setw(15)<<epsr[2][i]<<std::setprecision(2)<<std::setw(15)<<epsr[3][i]<<std::setprecision(2)<<std::setw(15)<<epsr[4][i]<<std::setprecision(2)<<std::setw(15)<<epsr[5][i]<<std::setprecision(2)<<std::setw(15)<<epsr[6][i]<<std::setprecision(2)<<std::setw(15)<<epsr[7][i]<<std::setprecision(2)<<std::setw(15)<<epsr[8][i]<<std::endl;
		}

		GlobalV::ofs_running<<std::endl;
		GlobalV::ofs_running<<" The imag part of the macroscopic dielectric constant:"<<std::endl;
		GlobalV::ofs_running<<std::setw(15)<<"omega"<<std::setw(15)<<"XX"<<std::setw(15)<<"XY"<<std::setw(15)<<"XZ"<<std::setw(15)<<"YX"<<std::setw(15)<<"YY"<<std::setw(15)<<"YZ"<<std::setw(15)<<"ZX"<<std::setw(15)<<"ZY"<<std::setw(15)<<"ZZ"<<std::endl;
		for(int i=0; i<nomega; i++)
		{
			GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<(i*domega)<<std::setprecision(2)<<std::setw(15)<<epsi[0][i]<<std::setprecision(2)<<std::setw(15)<<epsi[1][i]<<std::setprecision(2)<<std::setw(15)<<epsi[2][i]<<std::setprecision(2)<<std::setw(15)<<epsi[3][i]<<std::setprecision(2)<<std::setw(15)<<epsi[4][i]<<std::setprecision(2)<<std::setw(15)<<epsi[5][i]<<std::setprecision(2)<<std::setw(15)<<epsi[6][i]<<std::setprecision(2)<<std::setw(15)<<epsi[7][i]<<std::setprecision(2)<<std::setw(15)<<epsi[8][i]<<std::endl;
		}

		/*for(int iw=0; iw<nomega; iw++)
		{
			std::cout <<"epsi[0]["<<iw<<"] = "<< epsi[0][iw]<<std::endl;
			std::cout <<"epsr[0]["<<iw<<"] = "<< epsr[0][iw]<<std::endl;
		}*/

		GlobalV::ofs_running<<" Macroscopic dielectric constant matrix :"<<std::endl;
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<epsr[0][0]<<std::setprecision(2)<<std::setw(15)<<epsr[1][0]<<std::setprecision(2)<<std::setw(15)<<epsr[2][0]<<std::endl;
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<epsr[3][0]<<std::setprecision(2)<<std::setw(15)<<epsr[4][0]<<std::setprecision(2)<<std::setw(15)<<epsr[5][0]<<std::endl;
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<epsr[6][0]<<std::setprecision(2)<<std::setw(15)<<epsr[7][0]<<std::setprecision(2)<<std::setw(15)<<epsr[8][0]<<std::endl;

		//std::cout <<"Macroscopic dielectric constant matrix :"<<std::endl;
		//std::cout << epsr[0][0] <<"  "<<epsr[1][0] <<"  "<<epsr[2][0]<<"  "<<std::endl;
		//std::cout << epsr[3][0] <<"  "<<epsr[4][0] <<"  "<<epsr[5][0]<<"  "<<std::endl;
		//std::cout << epsr[6][0] <<"  "<<epsr[7][0] <<"  "<<epsr[8][0]<<"  "<<std::endl;
	}

	Delete();

	return;
}

void Epsilon0_pwscf:: Init()
{

	dipole_aux = new std::complex<double>**[3];
	for(int i=0; i<3; i++)
	{
		dipole_aux[i] = new std::complex<double>*[GlobalV::NBANDS];
		for(int j=0; j<GlobalV::NBANDS; j++)
		{
			dipole_aux[i][j] = new std::complex<double>[GlobalV::NBANDS];
		}
	}

	dipole = new std::complex<double>**[9];
	for(int i=0; i<9; i++)
	{
		dipole[i] = new std::complex<double>*[GlobalV::NBANDS];
		for(int j=0; j<GlobalV::NBANDS; j++)
		{
			dipole[i][j] = new std::complex<double>[GlobalV::NBANDS];
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
			for(int j=0; j<GlobalV::NBANDS; j++)
			{
				delete[] dipole_aux[i][j];
			}
			delete[] dipole_aux[i];
		}
		delete[] dipole_aux;

		for(int i=0; i<9; i++)
		{
			for(int j=0; j<GlobalV::NBANDS; j++)
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
	if(GlobalV::NSPIN == 1)
		return (GlobalC::wf.wg(ik,ib) * 2.0/ GlobalC::kv.wk[ik] );
	else if(GlobalV::NSPIN == 2)
		return (GlobalC::wf.wg(ik,ib) * 1.0/ GlobalC::kv.wk[ik] );
	else				// Peize Lin add 2019-05-01
		throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
} 
void Epsilon0_pwscf:: Cal_dipole(int ik)
{

	std::complex<double> dipole_aux_core[3][GlobalV::NBANDS][GlobalV::NBANDS];

	for(int i=0; i<3; i++)
		for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				dipole_aux_core[i][ib1][ib2] = std::complex<double>(0.0,0.0);
			}


	for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
	{
		if( focc(ib2,ik) < 2.0)
		{
			for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
			{
				if(ib1 == ib2)
					continue;
				if( focc(ib1,ik) > 0.0001)
				{
					for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
					{
						dipole_aux_core[0][ib1][ib2] += conj(GlobalC::wf.evc[ik](ib1, ig)) * GlobalC::pw.get_G_cartesian_projection(ig, 0) * (ModuleBase::TWO_PI / GlobalC::ucell.lat0) * GlobalC::wf.evc[ik](ib2, ig);
						dipole_aux_core[1][ib1][ib2] += conj(GlobalC::wf.evc[ik](ib1, ig)) * GlobalC::pw.get_G_cartesian_projection(ig, 1) * (ModuleBase::TWO_PI / GlobalC::ucell.lat0) * GlobalC::wf.evc[ik](ib2,ig);
						dipole_aux_core[2][ib1][ib2] += conj(GlobalC::wf.evc[ik](ib1, ig)) * GlobalC::pw.get_G_cartesian_projection(ig, 2) * (ModuleBase::TWO_PI / GlobalC::ucell.lat0) * GlobalC::wf.evc[ik](ib2,ig);
					}
				}
			}

		}
	}

	if(metalcalc)
	{
		for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
		{
			for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
			{
				dipole_aux_core[0][ib1][ib1] += conj(GlobalC::wf.evc[ik](ib1, ig)) * GlobalC::pw.get_GPlusK_cartesian_projection(ik, ig, 0) * (ModuleBase::TWO_PI / GlobalC::ucell.lat0) * GlobalC::wf.evc[ik](ib1,ig);
				dipole_aux_core[1][ib1][ib1] += conj(GlobalC::wf.evc[ik](ib1, ig)) * GlobalC::pw.get_GPlusK_cartesian_projection(ik, ig, 1) * (ModuleBase::TWO_PI / GlobalC::ucell.lat0) * GlobalC::wf.evc[ik](ib1, ig);
				dipole_aux_core[2][ib1][ib1] += conj(GlobalC::wf.evc[ik](ib1, ig)) * GlobalC::pw.get_GPlusK_cartesian_projection(ik, ig, 2) * (ModuleBase::TWO_PI / GlobalC::ucell.lat0) * GlobalC::wf.evc[ik](ib1, ig);
			}
		}
	}

	double dipole_aux_core_R[3][GlobalV::NBANDS][GlobalV::NBANDS];
	double dipole_aux_core_I[3][GlobalV::NBANDS][GlobalV::NBANDS];
	double dipole_aux_R[3][GlobalV::NBANDS][GlobalV::NBANDS];
	double dipole_aux_I[3][GlobalV::NBANDS][GlobalV::NBANDS];

	for(int i=0; i<3; i++)
		for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				dipole_aux_core_R[i][ib1][ib2] = dipole_aux_core[i][ib1][ib2].real();
				dipole_aux_core_I[i][ib1][ib2] = dipole_aux_core[i][ib1][ib2].imag();
			}

#ifdef __MPI
	MPI_Allreduce(dipole_aux_core_R,dipole_aux_R,3*GlobalV::NBANDS*GlobalV::NBANDS,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	MPI_Allreduce(dipole_aux_core_I,dipole_aux_I,3*GlobalV::NBANDS*GlobalV::NBANDS,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif

	for(int i=0; i<3; i++)
		for(int ib1=0; ib1<GlobalV::NBANDS; ib1++)
			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				dipole_aux[i][ib1][ib2] = std::complex<double>( dipole_aux_R[i][ib1][ib2], dipole_aux_I[i][ib1][ib2]);
			}


	return;
}
