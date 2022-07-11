#include "./stress_func.h"
#include "../module_xc/xc_functional.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate the mGGA stress correction in PW and LCAO
void Stress_Func::stress_mgga(ModuleBase::matrix& sigma, const psi::Psi<complex<double>>* psi_in) 
{
	ModuleBase::timer::tick("Stress_Func","stress_mgga");

	if (GlobalV::NSPIN==4) ModuleBase::WARNING_QUIT("stress_mgga","noncollinear stress + mGGA not implemented");

	int current_spin = 0;
	
	std::complex<double>** gradwfc;
	std::complex<double>* psi;

	double*** crosstaus;

	int ipol2xy[3][3];
	double sigma_mgga[3][3];

	gradwfc = new std::complex<double>*[GlobalC::wfcpw->nrxx];
	crosstaus = new double**[GlobalC::wfcpw->nrxx];
	
	for(int ir = 0;ir<GlobalC::wfcpw->nrxx;ir++)
	{
		crosstaus[ir] = new double*[6];
		gradwfc[ir] = new std::complex<double>[3];
		ModuleBase::GlobalFunc::ZEROS(gradwfc[ir],3);
		for(int j = 0;j<6;j++)
		{
			crosstaus[ir][j] = new double [GlobalV::NSPIN];
			ModuleBase::GlobalFunc::ZEROS(crosstaus[ir][j],GlobalV::NSPIN);
		}
	}

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		if(GlobalV::NSPIN==2) current_spin = GlobalC::kv.isk[ik];
		const int npw = GlobalC::kv.ngk[ik]; 	
		psi = new complex<double>[npw];

		for (int ibnd = 0; ibnd < GlobalV::NBANDS; ibnd++)
		{
			const double w1 = GlobalC::wf.wg(ik, ibnd) / GlobalC::ucell.omega;
			const std::complex<double>* ppsi=nullptr;
			if(psi_in!=nullptr)
			{
				ppsi = &(psi_in[0](ik, ibnd, 0));
			}
			else
			{
				ppsi = &(GlobalC::wf.evc[ik](ibnd, 0));
			}
			for(int ig = 0; ig<npw; ig++)
			{
				psi[ig] = ppsi[ig];
			}
			XC_Functional::grad_wfc(psi, ik, gradwfc, GlobalC::wfcpw);

			int ipol = 0;
			for (int ix = 0; ix < 3; ix++)
			{
				for (int iy = 0; iy < ix+1; iy++)
				{
					ipol2xy[ix][iy]=ipol;
					ipol2xy[iy][ix]=ipol;
					for(int ir = 0;ir<GlobalC::wfcpw->nrxx;ir++)
					{
						crosstaus[ir][ipol][current_spin] += 2.0 * w1 * (gradwfc[ir][ix].real() * gradwfc[ir][iy].real() + gradwfc[ir][ix].imag() * gradwfc[ir][iy].imag());
					}
					ipol+=1;
				}
			}
		}//band loop
		delete[] psi;
	}//k loop

	//if we are using kpools, then there should be a 
	//reduction of crosstaus w.r.t. kpools here.
	//will check later

	for(int ir = 0;ir<GlobalC::wfcpw->nrxx;ir++)
	{
		delete[] gradwfc[ir];
	}
	delete[] gradwfc;

	for(int is = 0; is < GlobalV::NSPIN; is++)
	{
		for (int ix = 0; ix < 3; ix++)
		{
			for (int iy = 0; iy < 3; iy++)
			{
				double delta= 0.0;
				if(ix==iy) delta=1.0;
				sigma_mgga[ix][iy] = 0.0;
				for(int ir = 0;ir<GlobalC::wfcpw->nrxx;ir++)
				{
					double x = GlobalC::pot.vofk(is,ir) * (GlobalC::CHR.kin_r[is][ir] * delta + crosstaus[ir][ipol2xy[ix][iy]][is]);
					sigma_mgga[ix][iy] += x;
				}
			}
		}
	}
	
	for(int ir = 0;ir<GlobalC::wfcpw->nrxx;ir++)
	{
		for(int j = 0;j<6;j++)
		{
			delete[] crosstaus[ir][j];
		}
		delete[] crosstaus[ir];
	}
	delete[] crosstaus;

#ifdef __MPI
	for(int l = 0;l<3;l++)
	{
		for(int m = 0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_all( sigma_mgga[l][m] );
		}
	}
#endif	
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			sigma(i,j) += sigma_mgga[i][j] / GlobalC::wfcpw->nxyz;
		}
	}
	ModuleBase::timer::tick("Stress_Func","stress_mgga");
	return;
}
