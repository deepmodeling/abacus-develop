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
#include "epsilon0_vasp.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace GlobalC
{
Epsilon0_vasp epsilon0_vasp;
}

Epsilon0_vasp::Epsilon0_vasp()
               :init_finish(false)
{
	epsilon = false;
	domega = 0.01;
	nomega = 300;
	eta = 0.01;
}

Epsilon0_vasp::~Epsilon0_vasp()
{
}

void Epsilon0_vasp::cal_epsilon0(const elecstate::ElecState* pelec)
{
	std::cout << "domega = "<<domega<<std::endl;
	std::cout << "nomega = "<<nomega<<std::endl;
	std::cout << "eta = "<<eta<<std::endl;
	
	double occupied_bands = static_cast<double>(GlobalV::nelec/ModuleBase::DEGSPIN);
	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0;
	}

	if(GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
	{
		oband = int(occupied_bands + 0.5);
	}
	else
	{
		oband = 2 * int(occupied_bands + 0.5);
	}
	uband = GlobalV::NBANDS - oband;

	std::cout << "oband = "<<oband<<std::endl;
	std::cout << "uband = "<<uband<<std::endl;

	if( !init_finish )
	{
		Init();
		init_finish = true;
	}

	//GlobalC::ppcell.init_vnl_alpha();

	Cal_epsilon0s(pelec);

	Cal_T();
	Cal_epsilon0();

	GlobalV::ofs_running<<std::endl;
	GlobalV::ofs_running<<" The real part of the macroscopic dielectric constant:"<<std::endl;
	GlobalV::ofs_running<<std::setw(15)<<"omega"<<std::setw(15)<<"XX"<<std::setw(15)<<"XY"<<std::setw(15)<<"XZ"<<std::setw(15)<<"YX"<<std::setw(15)<<"YY"<<std::setw(15)<<"YZ"<<std::setw(15)<<"ZX"<<std::setw(15)<<"ZY"<<std::setw(15)<<"ZZ"<<std::endl;
	for(int i=0; i<nomega; i++)
	{
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<(i*domega)<<std::setprecision(2)<<std::setw(15)<<eps0[0][i].real()+1.0<<std::setprecision(2)<<std::setw(15)<<eps0[1][i].real()<<std::setprecision(2)<<std::setw(15)<<eps0[2][i].real()<<std::setprecision(2)<<std::setw(15)<<eps0[3][i].real()<<std::setprecision(2)<<std::setw(15)<<eps0[4][i].real()+1.0<<std::setprecision(2)<<std::setw(15)<<eps0[5][i].real()<<std::setprecision(2)<<std::setw(15)<<eps0[6][i].real()<<std::setprecision(2)<<std::setw(15)<<eps0[7][i].real()<<std::setprecision(2)<<std::setw(15)<<eps0[8][i].real()+1.0<<std::endl;
	}

	GlobalV::ofs_running<<std::endl;
	GlobalV::ofs_running<<" The imag part of the macroscopic dielectric constant:"<<std::endl;
	GlobalV::ofs_running<<std::setw(15)<<"omega"<<std::setw(15)<<"XX"<<std::setw(15)<<"XY"<<std::setw(15)<<"XZ"<<std::setw(15)<<"YX"<<std::setw(15)<<"YY"<<std::setw(15)<<"YZ"<<std::setw(15)<<"ZX"<<std::setw(15)<<"ZY"<<std::setw(15)<<"ZZ"<<std::endl;
	for(int i=0; i<nomega; i++)
	{
		GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<(i*domega)<<std::setprecision(2)<<std::setw(15)<<-eps0[0][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[1][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[2][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[3][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[4][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[5][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[6][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[7][i].imag()<<std::setprecision(2)<<std::setw(15)<<-eps0[8][i].imag()<<std::endl;
	}

	GlobalV::ofs_running<<" Macroscopic dielectric constant matrix :"<<std::endl;
	GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<eps0[0][0].real()+1.0<<std::setprecision(2)<<std::setw(15)<<eps0[1][0].real()<<std::setprecision(2)<<std::setw(15)<<eps0[2][0].real()<<std::endl;
	GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<eps0[3][0].real()<<std::setprecision(2)<<std::setw(15)<<eps0[4][0].real()+1.0<<std::setprecision(2)<<std::setw(15)<<eps0[5][0].real()<<std::endl;
	GlobalV::ofs_running<<std::setprecision(2)<<std::setw(15)<<eps0[6][0].real()<<std::setprecision(2)<<std::setw(15)<<eps0[7][0].real()<<std::setprecision(2)<<std::setw(15)<<eps0[8][0].real()+1.0<<std::endl;

	//std::cout <<"Macroscopic dielectric constant matrix :"<<std::endl;
	//std::cout << eps0[0][0].real()+1.0 <<"  "<<eps0[1][0].real() <<"  "<<eps0[2][0].real() <<"  "<<std::endl;
	//std::cout << eps0[3][0].real() <<"  "<<eps0[4][0].real()+1.0 <<"  "<<eps0[5][0].real() <<"  "<<std::endl;
	//std::cout << eps0[6][0].real() <<"  "<<eps0[7][0].real() <<"  "<<eps0[8][0].real()+1.0 <<"  "<<std::endl;

	Delete();

	return;
}

void Epsilon0_vasp:: Init()
{
	b = new std::complex<double> **[oband];
	for(int ib1=0; ib1<oband; ib1++)
	{
		b[ib1] = new std::complex<double> *[uband];
		for(int ib2=0; ib2<uband; ib2++)
		{
			b[ib1][ib2] = new std::complex<double>[3];
		}
	}

	psi = new std::complex<double> *[GlobalV::NBANDS];
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		psi[ib] = new std::complex<double>[GlobalC::wfcpw->nrxx];
	}

	psi_nabla =new std::complex<double> **[GlobalV::NBANDS];
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		psi_nabla[ib] = new std::complex<double> *[GlobalC::wfcpw->nrxx];
		for(int ir=0; ir<GlobalC::wfcpw->nrxx; ir++)
		{
			psi_nabla[ib][ir] = new std::complex<double>[3];
		}
	}

	psi_nu = new std::complex<double> **[GlobalV::NBANDS];
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		psi_nu[ib] = new std::complex<double> *[GlobalC::ppcell.nkb];
		for(int u=0; u<GlobalC::ppcell.nkb; u++)
		{
			psi_nu[ib][u] = new std::complex<double>[4];
		}
	}

	eps0s = new std::complex<double> *[9];
	for(int i=0; i<9; i++)
	{
		eps0s[i] = new std::complex<double>[nomega];
	}

	T = new std::complex<double> *[nomega];
	for(int i=0; i<nomega; i++)
	{
		T[i] = new std::complex<double>[nomega];
	}

	eps0 = new std::complex<double> *[9];
	for(int i=0; i<9; i++)
	{
		eps0[i] = new std::complex<double>[nomega];
	}

	return;
}

void Epsilon0_vasp:: Delete()
{
	for(int ib1=0; ib1<oband; ib1++)
	{
		for(int ib2=0; ib2<uband; ib2++)
		{
			delete[] b[ib1][ib2];
		}
		delete[] b[ib1];
	}
	delete[] b;

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		delete[] psi[ib];
	}
	delete[] psi;

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		for(int ir=0; ir<GlobalC::wfcpw->nrxx; ir++)
		{
			delete[] psi_nabla[ib][ir];
		}
		delete[] psi_nabla[ib];
	}
	delete[] psi_nabla;

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		for(int u=0; u<GlobalC::ppcell.nkb; u++)
		{
			delete[] psi_nu[ib][u];
		}
		delete[] psi_nu[ib];
	}
	delete[] psi_nu;

	for(int i=0; i<9; i++)
	{
		delete eps0s[i];
	}
	delete[] eps0s;

	for(int i=0; i<nomega; i++)
	{
		delete[] T[i];
	}
	delete[] T;

	for(int i=0; i<9; i++)
	{
		delete eps0[i];
	}
	delete[] eps0;

	return;
}

void Epsilon0_vasp:: Cal_psi(int ik)      // pengfei Li 2018-11-13
{
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		GlobalC::wfcpw->recip2real(&GlobalC::wf.evc[ik](ib, 0), psi[ib], ik);
	}
	return;
}

void Epsilon0_vasp:: Cal_psi_nabla(int ik)      // pengfei Li 2018-11-13
{
	std::complex<double>* porter = new std::complex<double>[GlobalC::wfcpw->nmaxgr];
	for(int i = 0 ; i < 3 ; ++i)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ++ib)
		{

			for(int ig = 0; ig < GlobalC::kv.ngk[ik] ; ig++)
			{
				porter[ig] = GlobalC::wf.evc[ik](ib,ig) * (GlobalC::wfcpw->getgpluskcar(ik, ig)[i] * (ModuleBase::TWO_PI/GlobalC::ucell.lat0));
			}

			GlobalC::wfcpw->recip2real(porter, porter, ik);

			for(int ir=0; ir<GlobalC::wfcpw->nrxx; ir++)
			{
				psi_nabla[ib][ir][i] = porter[ir];
			}
		}
	}
		
	delete[] porter;

	return;
}

void Epsilon0_vasp:: Cal_b(int ik)
{
	std::complex<double> b_core[oband][uband][3];

	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int i=0; i<3; i++)
			{
				b_core[ib1][ib2][i] = std::complex<double>(0.0,0.0);
			}

	//Cal_psi(ik);
	//Cal_psi_nabla(ik);

	if(GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
	{
		for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				for(int ig =0; ig< GlobalC::kv.ngk[ik]; ig++)
				{
					b_core[ib1][ib2][0] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig)[0] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
					b_core[ib1][ib2][1] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig)[1] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
					b_core[ib1][ib2][2] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig)[2] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);

				}
			}
	}
	else if(GlobalV::NSPIN == 4)
	{
		for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				for(int ig =0; ig< GlobalC::kv.ngk[ik]; ig++)
				{
					b_core[ib1][ib2][0] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig)[0] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
					b_core[ib1][ib2][1] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig)[1] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
					b_core[ib1][ib2][2] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig)[2] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
				}
			}

		for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				for(int ig =GlobalC::wf.npwx; ig< GlobalC::wf.npwx + GlobalC::kv.ngk[ik]; ig++)
				{
					b_core[ib1][ib2][0] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig-GlobalC::wf.npwx)[0] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
					b_core[ib1][ib2][1] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig-GlobalC::wf.npwx)[1] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
					b_core[ib1][ib2][2] += conj(GlobalC::wf.evc[ik](ib1, ig)) * (GlobalC::wfcpw->getgpluskcar(ik,ig-GlobalC::wf.npwx)[2] * (ModuleBase::TWO_PI / GlobalC::ucell.lat0)) * GlobalC::wf.evc[ik](oband + ib2, ig);
				}
			}
	}

	double b_core_R[oband][uband][3];
	double b_core_I[oband][uband][3];
	double b_R[oband][uband][3];
	double b_I[oband][uband][3];

	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int i=0; i<3; i++)
			{
				b_core_R[ib1][ib2][i] = b_core[ib1][ib2][i].real();
				b_core_I[ib1][ib2][i] = b_core[ib1][ib2][i].imag();
			}

#ifdef __MPI
	MPI_Allreduce(b_core_R,b_R,3*oband*uband,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	MPI_Allreduce(b_core_I,b_I,3*oband*uband,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif

	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int i=0; i<3; i++)
			{
				b[ib1][ib2][i] = std::complex<double>( b_R[ib1][ib2][i], b_I[ib1][ib2][i]);
			}

	//Cal_psi_nu(ik);                       // pengfei Li 2018-11-13

	/*for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
		{
			int count = 0;
			for(int iat=0; iat<GlobalC::ucell.nat; iat++)
			{
				int it =  GlobalC::ucell.iat2it[iat];
				int nht = GlobalC::ucell.atoms[it].nh;
				int count1 = count;
				for(int ih=0; ih<nht; ih++)
				{
					int ikb = count1 + ih;
					for(int jh=0; jh<nht; jh++)
					{
						int jkb = count1 + jh;
						b[ib1][ib2][0] += IMAG_UNIT * GlobalC::ppcell.deeq(0, iat, ih, jh) *
						(psi_nu[ib1][ikb][3] * conj(psi_nu[oband+ib2][jkb][0]) -  psi_nu[ib1][ikb][0] * conj(psi_nu[oband+ib2][jkb][3]) );
						b[ib1][ib2][1] += IMAG_UNIT * GlobalC::ppcell.deeq(0, iat, ih, jh) *
						(psi_nu[ib1][ikb][3] * conj(psi_nu[oband+ib2][jkb][1]) -  psi_nu[ib1][ikb][1] * conj(psi_nu[oband+ib2][jkb][3]) );
						b[ib1][ib2][2] += IMAG_UNIT * GlobalC::ppcell.deeq(0, iat, ih, jh) *
						(psi_nu[ib1][ikb][3] * conj(psi_nu[oband+ib2][jkb][2]) -  psi_nu[ib1][ikb][2] * conj(psi_nu[oband+ib2][jkb][3]) );

					}
					count++;
				}
			}
		}*/

	return;
}

#ifdef __LCAO
void Epsilon0_vasp:: Cal_psi_nu(int ik)
{
	GlobalC::ppcell.getvnl_alpha(ik);

	std::complex<double> psi_nu_core[GlobalV::NBANDS][GlobalC::ppcell.nkb][4];

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
		for(int ikb=0; ikb<GlobalC::ppcell.nkb; ikb++)
			for(int i=0; i<4; i++)
			{
				psi_nu_core[ib][ikb][i] = std::complex<double>(0.0,0.0);
			}

	/*for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		for(int u=0; u<GlobalC::ppcell.nkb; u++)
		{
			for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
			{
				std::cout<<"ib = "<<ib<<" u = "<<u<<" ig = "<<ig<<" vkb_alpha = "<<GlobalC::ppcell.vkb_alpha[0][u][ig]<<std::endl;
			}
		}
	}*/

	/*for(int u=0; u<GlobalC::ppcell.nkb; u++)
	{
		for(int ig=0; ig<GlobalC::kv.ngk[ik];ig++)
		{
			std::cout<<
		}
	}*/

	for(int u=0; u<GlobalC::ppcell.nkb; u++)
		for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
		{
			std::cout<<"ik = "<<ik<<" u = "<<u<<" ig = "<<ig<<"  GlobalC::ppcell.vkb["<<u<<"]["<<ig<<"] = "<<GlobalC::ppcell.vkb(u,ig)<<std::endl;
		}

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		for(int u=0; u<GlobalC::ppcell.nkb; u++)
		{
			ModuleBase::GlobalFunc::ZEROS(psi_nu[ib][u], 4);
			for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
			{
				//std::cout<<"ib = "<<ib<<" u = "<<u<<" ig = "<<ig<<" vkb_alpha = "<<GlobalC::ppcell.vkb_alpha[0][u][ig]<<std::endl;
				//std::cout<<"psi_nu = "<<psi_nu[ib][u][0]<<std::endl;
				psi_nu_core[ib][u][0] += conj(GlobalC::wf.evc[ik](ib,ig)) * GlobalC::ppcell.vkb_alpha[0][u][ig];
				psi_nu_core[ib][u][1] += conj(GlobalC::wf.evc[ik](ib,ig)) * GlobalC::ppcell.vkb_alpha[1][u][ig];
				psi_nu_core[ib][u][2] += conj(GlobalC::wf.evc[ik](ib,ig)) * GlobalC::ppcell.vkb_alpha[2][u][ig];
				psi_nu_core[ib][u][3] += conj(GlobalC::wf.evc[ik](ib,ig)) * GlobalC::ppcell.vkb(u,ig);
			}
		}
	}

	double psi_nu_core_R[GlobalV::NBANDS][GlobalC::ppcell.nkb][4];
	double psi_nu_core_I[GlobalV::NBANDS][GlobalC::ppcell.nkb][4];
	double psi_nu_R[GlobalV::NBANDS][GlobalC::ppcell.nkb][4];
	double psi_nu_I[GlobalV::NBANDS][GlobalC::ppcell.nkb][4];

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
		for(int ikb=0; ikb<GlobalC::ppcell.nkb; ikb++)
			for(int i=0; i<4; i++)
			{
				psi_nu_core_R[ib][ikb][i] = psi_nu_core[ib][ikb][i].real();
				psi_nu_core_I[ib][ikb][i] = psi_nu_core[ib][ikb][i].imag();
			}

#ifdef __MPI
	MPI_Allreduce(psi_nu_core_R,psi_nu_R,4*GlobalV::NBANDS*GlobalC::ppcell.nkb,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	MPI_Allreduce(psi_nu_core_I,psi_nu_I,4*GlobalV::NBANDS*GlobalC::ppcell.nkb,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
		for(int ikb=0; ikb<GlobalC::ppcell.nkb; ikb++)
			for(int i=0; i<4; i++)
			{
				psi_nu[ib][ikb][i] = std::complex<double>( psi_nu_R[ib][ikb][i], psi_nu_I[ib][ikb][i]);
			}

	return;
}
#endif

void Epsilon0_vasp:: Cal_epsilon0s(const elecstate::ElecState* pelec)
{
	for(int i=0; i<9; i++)
		for(int j=0; j<nomega; j++)
		{
			eps0s[i][j] = std::complex<double>(0.0, 0.0);
		}

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		Cal_b(ik);
		/*for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				std::cout<<"ik = "<<ik<<" b["<<ib1<<"]["<<ib2<<"][0]"<<b[ib1][ib2][0]<<std::endl;
			}*/

		for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				double delta_e = pelec->ekb(ik, oband+ib2) - pelec->ekb(ik, ib1);
				if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
				{
					int n = int(delta_e/domega);
					double e1 = double(n) * domega;
					double e2 = double(n+1) * domega;
					std::complex<double> weight1 = std::complex<double>(32 * (pelec->wg(ik,ib1) - pelec->wg(ik,oband+ib2)) * (e2 - delta_e)/domega/GlobalC::ucell.omega * ModuleBase::PI * ModuleBase::PI /delta_e/delta_e/domega, 0.0);
					std::complex<double> weight2 = std::complex<double>(32 * (pelec->wg(ik,ib1) - pelec->wg(ik,oband+ib2)) * (delta_e - e1)/domega/GlobalC::ucell.omega * ModuleBase::PI * ModuleBase::PI /delta_e/delta_e/domega, 0.0);
					
					eps0s[0][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][0]);
					eps0s[1][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][1]);
					eps0s[2][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][2]);
					eps0s[3][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][0]);
					eps0s[4][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][1]);
					eps0s[5][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][2]);
					eps0s[6][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][0]);
					eps0s[7][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][1]);
					eps0s[8][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][2]);

					eps0s[0][n+1] += weight2 * b[ib1][ib2][0] * conj(b[ib1][ib2][0]);
					eps0s[1][n+1] += weight2 * b[ib1][ib2][0] * conj(b[ib1][ib2][1]);
					eps0s[2][n+1] += weight2 * b[ib1][ib2][0] * conj(b[ib1][ib2][2]);
					eps0s[3][n+1] += weight2 * b[ib1][ib2][1] * conj(b[ib1][ib2][0]);
					eps0s[4][n+1] += weight2 * b[ib1][ib2][1] * conj(b[ib1][ib2][1]);
					eps0s[5][n+1] += weight2 * b[ib1][ib2][1] * conj(b[ib1][ib2][2]);
					eps0s[6][n+1] += weight2 * b[ib1][ib2][2] * conj(b[ib1][ib2][0]);
					eps0s[7][n+1] += weight2 * b[ib1][ib2][2] * conj(b[ib1][ib2][1]);
					eps0s[8][n+1] += weight2 * b[ib1][ib2][2] * conj(b[ib1][ib2][2]);

				}
				else if((delta_e > ((nomega-1) * domega) || delta_e == ((nomega-1) * domega)) && delta_e < (nomega * domega))
				{
					int n = int(delta_e/domega);
					double e1 = double(n) * domega;
					double e2 = double(n+1) * domega;					
					std::complex<double> weight1 = std::complex<double>(32 * (pelec->wg(ik,ib1) - pelec->wg(ik,oband+ib2)) * (e2 - delta_e)/domega/GlobalC::ucell.omega * ModuleBase::PI * ModuleBase::PI /delta_e/delta_e/domega, 0.0);

					eps0s[0][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][0]);
					eps0s[1][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][1]);
					eps0s[2][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][2]);
					eps0s[3][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][0]);
					eps0s[4][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][1]);
					eps0s[5][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][2]);
					eps0s[6][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][0]);
					eps0s[7][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][1]);
					eps0s[8][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][2]);

				}

			}
	}

	return;
}

void Epsilon0_vasp:: Cal_T()
{
	std::complex<double> M1, M2;
	for(int n1=0; n1<nomega; n1++)
		for(int n=0; n<nomega; n++)
		{
			double n_e = double(n) * domega;
			double n1_e = double(n1) * domega;

			M1 = std::complex<double>(n1_e - n_e, -eta);
			M2 = std::complex<double>(-n1_e - n_e, eta);
			
			T[n][n1] = (1.0/M1 + 1.0/M2) * (-0.5) * 2.0/ModuleBase::PI * domega;
		}

	return;
}

void Epsilon0_vasp:: Cal_epsilon0()
{
	CMatrixMul(9,nomega,nomega,eps0s,T,eps0);

	return;
}

void Epsilon0_vasp:: CMatrixMul(int m, int n, int l, std::complex<double>** A, std::complex<double>** B, std::complex<double>** C)
{
	for(int i=0; i<m; i++)
		for(int j=0; j<l; j++)
		{
			C[i][j] = std::complex<double>(0.0, 0.0);
			for(int k=0; k<n; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}

	return;
}
