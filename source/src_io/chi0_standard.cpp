//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-30
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating some optical properties using the linear response
//theory .
//-------------------------------------------------------------------------------

#include "../src_pw/global.h"
#include "chi0_standard.h"
#include "../src_pw/hamilt.h"
#include "../src_lcao/wavefunc_in_pw.h"
#include "../src_pw/klist.h"
#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

namespace GlobalC
{
Chi0_standard chi0_standard;
}

Chi0_standard:: Chi0_standard()
                :init_finish(false)
{
      epsilon = false;
	  system = "bulk";
	  eta = 0.05;
	  domega = 0.01;
	  nomega = 300;
	  dim = 1;
	  oband = 1;
	  start_q = 1;
	  interval_q = 1;
	  nq = 1;
	  out_epsilon = true;
}

Chi0_standard:: ~Chi0_standard()
{
}

void Chi0_standard:: Chi(const elecstate::ElecState* pelec)
{
	ModuleBase::TITLE("Chi0_standard","Chi");

	//---------------------------------------
	//  the serial number of q
	//---------------------------------------
	bool exist_q = false;
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		double dx = fabs( GlobalC::kv.kvec_d[ik].x - q_start[0]);
		double dy = fabs( GlobalC::kv.kvec_d[ik].y - q_start[1]);
		double dz = fabs( GlobalC::kv.kvec_d[ik].z - q_start[2]);
		if( dx<0.0001 && dy<0.0001 && dz<0.0001)
		{
			start_q = ik; exist_q = true;
			break;
		}
	}

	if(!exist_q)
	{
		ModuleBase::WARNING_QUIT("chi0_hilbert","the chosen q is not included in the kmesh!!");
	}
	else
	{
		std::cout <<"start_q = "<<start_q<<std::endl;
	}


	int icount = 0;
	std::vector<int> qc(GlobalC::kv.nks);		// Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<double> ql(GlobalC::kv.nks);		// Peize Lin change ptr to std::vector at 2020.01.31
	int total_icount=0;
	int temp1; double temp2;

	if( direct[0]!=0 && direct[1]!=0 && direct[2]!=0)
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p0 = x/direct[0]; double p1 = y/direct[1]; double p2 = z/direct[2];
			if( p0>0.0001 && fabs(p0-p1)<0.0001 && fabs(p0-p2)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}
	}
	else if( direct[0]==0 && direct[1]!=0 && direct[2]!=0 )
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p1 = y/direct[1]; double p2 = z/direct[2];
			if( fabs(x)<0.0001 && p1>0.0001 && fabs(p1-p2)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p1; icount++;
			}
		}
	}
	else if( direct[0]!=0 && direct[1]==0 && direct[2]!=0 )
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p0 = x/direct[0]; double p2 = z/direct[2];
			if( fabs(y)<0.0001 && p0>0.0001 && fabs(p0-p2)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}
	}
	else if( direct[0]!=0 && direct[1]!=0 && direct[2]==0 )
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p0 = x/direct[0]; double p1 = y/direct[1];
			if( fabs(z)<0.0001 && p0>0.0001 && fabs(p0-p1)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}
	}
	else if( direct[0]==0 && direct[1]==0 && direct[2]!=0 )
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p2 = z/direct[2];
			if( fabs(x)<0.0001 && fabs(y)<0.0001 && p2 > 0.0001)
			{
				qc[icount] = ik; ql[icount] = p2; icount++;
			}
		}
	}
	else if( direct[0]==0 && direct[1]!=0 && direct[2]==0 )
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p1 = y/direct[1];
			if( fabs(x)<0.0001 && fabs(z)<0.0001 && p1 > 0.0001)
			{
				qc[icount] = ik; ql[icount] = p1; icount++;
			}
		}
	}
	else if( direct[0]!=0 && direct[1]==0 && direct[2]==0 )
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			double x = GlobalC::kv.kvec_d[ik].x - GlobalC::kv.kvec_d[start_q].x;
			double y = GlobalC::kv.kvec_d[ik].y - GlobalC::kv.kvec_d[start_q].y;
			double z = GlobalC::kv.kvec_d[ik].z - GlobalC::kv.kvec_d[start_q].z;

			double p0 = x/direct[0];
			if( fabs(y)<0.0001 && fabs(z)<0.0001 && p0 > 0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}
	}
	total_icount = icount;

	//if(total_icount == 0)
	//{
	//	ModuleBase::WARNING_QUIT("chi0_hilbert","Now the kmesh contains no kpoint along this direction!");
	//}
	if(total_icount < nq-1 )
	{
		ModuleBase::WARNING_QUIT("chi0_hilbert","Now the kmesh doesn't contain enough kpoints along this direction! please change the parameter nq smaller");
	}

	for(icount=0; icount<total_icount-1; icount++)
	{
		if(ql[icount] < ql[icount+1])
		{
			temp1 = qc[icount+1]; qc[icount+1] = qc[icount]; qc[icount] = temp1;
			temp2 = ql[icount+1]; ql[icount+1] = ql[icount]; ql[icount] = temp2;
		}
	}

	int temp_q = qc[total_icount-1];
	//std::cout <<"temp_q = "<<temp_q<<std::endl;

	interval_q = temp_q - start_q;
	//std::cout <<"interval_q = "<<interval_q<<std::endl;

	Parallel_G();

	if(system == "surface")
	{
		dim_para = parallel_g();
		std::cout << "dim_para = "<<dim_para<<std::endl;
	}

	//----------------------------------------------------------
	// calculate the number of occupied bands
	//----------------------------------------------------------
	int occ_bands = 0;
	bool occ_flag;
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		occ_flag = false;
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			if( pelec->wg(ik,ib)> 0.0001)
			{
				occ_flag = true;
				continue;
			}
		}

		if(occ_flag == true)
			occ_bands++;
	}
	std::cout <<"occ_bands = "<<occ_bands<<std::endl;
	oband = occ_bands;

	if(!init_finish)
	{
		Init();
		init_finish = true;
	}

	for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
	{
		double q = sqrt(((GlobalC::kv.kvec_c[iq])*(ModuleBase::TWO_PI/GlobalC::ucell.lat0)).norm2());
		double gather[nomega];

		int count =0;
		for(double omega=0.0; omega<(domega*nomega); omega=omega+domega)
		{
			Cal_chi0(iq,omega, pelec);
			std::cout<<"chi0 iq= "<<iq<<" omega= "<<omega<<"  "<<chi0[0][0].real()<<" "<<chi0[0][0].imag()<<std::endl;
			if(system == "surface")
			{
				chi0_para_g();
				std::cout<<"chi0_para iq= "<<iq<<" omega= "<<omega<<"  "<<chi0_para[0][0].real()<<" "<<chi0_para[0][0].imag()<<std::endl;
			}
			Cal_rpa(iq);
			Cal_chi();
			std::cout<<"chi iq= "<<iq<<" omega= "<<omega<<"  "<<chi[0][0].real()<<" "<<chi[0][0].imag()<<std::endl;
			std::cout<<"epsilon iq= "<<iq<<" omega= "<<omega<<"  "<<8*ModuleBase::PI/q/q*chi[0][0].real()<<" "<<8*ModuleBase::PI/q/q*chi[0][0].imag()<<std::endl;
			gather[count] = -8*ModuleBase::PI/q/q*chi[0][0].imag();		
			count++;
		}

		if(out_epsilon)
		{
			std::stringstream sseps;
			sseps << GlobalV::global_out_dir << "Imeps^-1_"<<iq<<".dat";
			std::ofstream ofseps(sseps.str().c_str());
			ofseps<<"Energy(Ry)"<<"   "<<"-Im{epsilon^-1}"<<std::endl;
			for(int i=0; i<nomega; i++)
			{
				ofseps <<i*domega<<"   "<<gather[i]<<std::endl;
			}
			ofseps.close();
		}
	}

	Delete();
	return;
}

void Chi0_standard::Parallel_G()
{
	ModuleBase::TITLE("Chi0_standard","Parallel_G");
	//----------------------------
	// init
	//----------------------------
	num_G_core = new int[GlobalV::DSIZE];
	num_G_dis = new int[GlobalV::DSIZE];
	G_r_core = new double[GlobalC::rhopw->npw];
	num_Gvector_core = new int[GlobalV::DSIZE];
	num_Gvector_dis = new int[GlobalV::DSIZE];
	G_r = new double[GlobalC::rhopw->npwtot];
	Gvec_core = new double[3*GlobalC::rhopw->npw];
	Gvec = new double[3*GlobalC::rhopw->npwtot];
	all_gcar = new ModuleBase::Vector3<double>[GlobalC::rhopw->npwtot];
	flag = new int[GlobalC::rhopw->npwtot];

	for(int i=0;i<GlobalC::rhopw->npwtot;i++)
	{
		flag[i] = i;
	}
	
	ModuleBase::GlobalFunc::ZEROS( num_G_dis, GlobalV::DSIZE);
	ModuleBase::GlobalFunc::ZEROS( num_G_core, GlobalV::DSIZE);
	ModuleBase::GlobalFunc::ZEROS( num_Gvector_dis, GlobalV::DSIZE);
	ModuleBase::GlobalFunc::ZEROS( num_Gvector_core, GlobalV::DSIZE);
	
#ifdef __MPI
	MPI_Allgather( &GlobalC::rhopw->npw, 1, MPI_INT, num_G_core, 1, MPI_INT, POOL_WORLD);
#endif

	memset(num_G_dis,0,GlobalV::DSIZE*sizeof(int));
	for(int i=0; i<GlobalV::DSIZE; i++)
	{
		for(int j=0; j<i; j++)
		{
			num_G_dis[i] += num_G_core[j];
		}
	}

	for(int i=0;i<GlobalV::DSIZE;i++)
	{
		num_Gvector_dis[i] = num_G_dis[i] * 3;
		num_Gvector_core[i] = num_G_core[i] * 3;
	}

	for(int g0=0;g0<GlobalC::rhopw->npw; g0++)
	{
		G_r_core[g0] = GlobalC::rhopw->gg[g0];
		Gvec_core[3 * g0] = GlobalC::rhopw->gcar[g0][0];
		Gvec_core[3 * g0 + 1] = GlobalC::rhopw->gcar[g0][1];
		Gvec_core[3 * g0 + 2] = GlobalC::rhopw->gcar[g0][2];
	}

#ifdef __MPI
	MPI_Allgatherv( G_r_core, GlobalC::rhopw->npw, MPI_DOUBLE, G_r, num_G_core, num_G_dis, MPI_DOUBLE, POOL_WORLD);
	MPI_Allgatherv( Gvec_core, 3*GlobalC::rhopw->npw, MPI_DOUBLE, Gvec, num_Gvector_core, num_Gvector_dis, MPI_DOUBLE, POOL_WORLD);
#endif

	double t1; int t2;
	for(int i=0;i<GlobalC::rhopw->npwtot;i++)
	{
		for(int j=0;j<GlobalC::rhopw->npwtot-i-1;j++)
		{
			if(G_r[j]>G_r[j+1])
			{
				t1=G_r[j]; G_r[j]=G_r[j+1]; G_r[j+1]=t1;
				t2=flag[j]; flag[j]=flag[j+1]; flag[j+1]=t2;
				t1=Gvec[3*j]; Gvec[3*j]=Gvec[3*j+3]; Gvec[3*j+3]=t1;
				t1=Gvec[3*j+1]; Gvec[3*j+1]=Gvec[3*j+4]; Gvec[3*j+4]=t1;
				t1=Gvec[3*j+2]; Gvec[3*j+2]=Gvec[3*j+5]; Gvec[3*j+5]=t1;
			}
		}
	}

	for(int i=0;i<GlobalC::rhopw->npwtot;i++)
	{
		all_gcar[i].x = Gvec[3*i]; all_gcar[i].y = Gvec[3*i+1]; all_gcar[i].z = Gvec[3*i+2];
		//std::cout<<"all_gcar["<<i<<"]= "<<all_gcar[i].x<<" "<<all_gcar[i].y<<" "<<all_gcar[i].z<<std::endl;
	}

	return;
}

void Chi0_standard:: Init()
{

	b_core = new std::complex<double>[GlobalC::rhopw->npw];

	b_summary = new std::complex<double>[GlobalC::rhopw->npwtot];

	b_order = new std::complex<double>[GlobalC::rhopw->npwtot];

	psi_r1 = new std::complex<double>*[GlobalV::NBANDS];
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		psi_r1[ib] = new std::complex<double>[GlobalC::wfcpw->nrxx];
	}

	psi_r2 = new std::complex<double>*[GlobalV::NBANDS];
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		psi_r2[ib] = new std::complex<double>[GlobalC::wfcpw->nrxx];
	}
	std::cout << "psi OK" <<std::endl;

	b = new std::complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		b[g0] = new std::complex<double>[oband*GlobalV::NBANDS];
	}
	std::cout << "b ok"<<std::endl;

	A = new std::complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		A[g0] = new std::complex<double>[oband*GlobalV::NBANDS];
	}
	std::cout << "A ok"<<std::endl;

	B = new std::complex<double>*[oband*GlobalV::NBANDS];
	for(int ib=0; ib<(oband*GlobalV::NBANDS); ib++)
	{
		B[ib] = new std::complex<double>[dim];
	}
	std::cout << "B ok" << std::endl;

	weight = new std::complex<double> [oband*GlobalV::NBANDS];
	std::cout << "weight OK"<<std::endl;

	chi0 = new std::complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		chi0[g0] = new std::complex<double>[dim];
	}
	std::cout << "chi0 OK" <<std::endl;

	chi0_para = new std::complex<double>*[dim_para];
	for(int g0=0; g0<dim_para; g0++)
	{
		chi0_para[g0] = new std::complex<double>[dim_para];
	}

	rpa = new std::complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		rpa[g0] = new std::complex<double>[dim];
	}
	std::cout << "rpa OK" <<std::endl;

	chi = new std::complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		chi[g0] = new std::complex<double>[dim];
	}
	std::cout << "chi OK" <<std::endl;

	std::cout << "All ok"<<std::endl;

	return;
}

void Chi0_standard::Delete()
{
	ModuleBase::TITLE("Chi0_standard","Delete");
	if(init_finish == true)				// Peize Lin change = to == at 2020.01.31
	{
		delete[] b_core;
		delete[] b_summary;
		delete[] b_order;

		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			delete[] psi_r1[ib];
		}
		delete[] psi_r1;

		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			delete[] psi_r2[ib];
		}
		delete[] psi_r2;

		for(int g0=0; g0<dim; g0++)
		{
			delete[] b[g0];
		}
		delete[] b;

		for(int g0=0; g0<dim; g0++)
		{
			delete[] A[g0];
		}
		delete[] A;

		for(int ib=0; ib<(oband*GlobalV::NBANDS); ib++)
		{
			delete[] B[ib];
		}
		delete[] B;

		delete[] weight;

		for(int g0=0; g0<dim; g0++)
		{
			delete[] chi0[g0];
		}
		delete[] chi0;

		for(int g0=0; g0<dim_para;g0++)
		{
			delete[] chi0_para[g0];
		}
		delete[] chi0_para;

		for(int g0=0; g0<dim; g0++)
		{
			delete[] rpa[g0];
		}
		delete[] rpa;

		for(int g0=0; g0<dim; g0++)
		{
			delete[] chi[g0];
		}
		delete[] chi;

	}
}

void Chi0_standard::Cal_Psi(int iq, std::complex<double> **psi_r)
{
	double phase_x, phase_xy, phase_xyz;
	std::complex<double> exp_tmp;
	std::complex<double>* porter = new std::complex<double> [GlobalC::wfcpw->nrxx];
	for(int ib = 0; ib < GlobalV::NBANDS; ib++)
	{

		GlobalC::wfcpw->recip2real(&GlobalC::wf.evc[iq](ib, 0), porter, iq);

		int ir=0;
		for(int ix=0; ix<GlobalC::wfcpw->nx; ix++)
		{
			phase_x = GlobalC::kv.kvec_d[iq].x*ix/GlobalC::wfcpw->nx;
			for(int iy=0; iy<GlobalC::wfcpw->ny; iy++)
			{
				phase_xy = phase_x + GlobalC::kv.kvec_d[iq].y*iy/GlobalC::wfcpw->ny;
				for(int iz=GlobalC::wfcpw->startz_current; iz<GlobalC::wfcpw->startz_current+GlobalC::wfcpw->nplane; iz++)
				{
					phase_xyz = (phase_xy + GlobalC::kv.kvec_d[iq].z*iz/GlobalC::wfcpw->nz) *ModuleBase::TWO_PI;
					exp_tmp = std::complex<double>( cos(phase_xyz), sin(phase_xyz) );
					psi_r[ib][ir] = porter[ir]*exp_tmp;
					ir++;
				}

			}
		}
	}
	delete[] porter;
	return;
}

void Chi0_standard::Cal_b(int iq, int ik, int iqk,  ModulePW::PW_Basis *rho_basis)
{
	ModuleBase::TITLE("Chi0_standard","Cal_b");
	ModuleBase::Vector3<double> qk;
	qk = GlobalC::kv.kvec_c[iq] + GlobalC::kv.kvec_c[ik];
	//std::cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<std::endl;
	double phase_x, phase_xy, phase_xyz;
	ModuleBase::Vector3<double> q = GlobalC::kv.kvec_d[iq];
	std::complex<double>* aux = new std::complex<double>[rho_basis->nmaxgr];

	Cal_Psi(ik, psi_r1);
	Cal_Psi(iqk, psi_r2);

	const int startz = rho_basis->startz_current; 
	const int nplane = rho_basis->nplane;

	for(int ib1=0; ib1< oband; ib1++)
	{
		for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
		{
			int ir=0;
			for(int ix=0; ix<rho_basis->nx; ix++)
			{
				phase_x = q.x*ix/rho_basis->nx;
				for(int iy=0; iy<rho_basis->ny; iy++)
				{
					phase_xy = phase_x + q.y*iy/rho_basis->ny;
					for(int iz=startz; iz<startz+nplane; iz++)
					{
						phase_xyz = (phase_xy + q.z*iz/rho_basis->nz) *ModuleBase::TWO_PI;
						std::complex<double> exp_tmp = std::complex<double>(cos(-phase_xyz), sin(-phase_xyz));
						aux[ir] = conj(psi_r1[ib1][ir]) * psi_r2[ib2][ir] *exp_tmp;
						ir++;
					}
				}
			}

			rho_basis->real2recip(aux,aux);

			for(int ig=0 ; ig < rho_basis->npw; ++ig)
			{
				b_core[ig] = aux[ ig ];
			}

#ifdef __MPI
			MPI_Allgatherv( b_core, rho_basis->npw, mpicomplex, b_summary, num_G_core, num_G_dis, mpicomplex, POOL_WORLD);
#endif

			for(int i=0;i<rho_basis->npwtot;i++)
			{
				b_order[i] = b_summary[flag[i]];
			}

			for(int g0=0; g0<dim; g0++)
			{
				b[g0][ib2+ib1*GlobalV::NBANDS] = b_order[g0];
			}

		}
	}

	delete[] aux;
	return;
}

void Chi0_standard:: Cal_weight(int iq, int ik, double omega, const elecstate::ElecState* pelec)
{
	int iqk = Cal_iq(ik, iq, GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
		{
			std::complex<double> factor = std::complex<double>( (omega + pelec->ekb(ik, ib1) - pelec->ekb(iqk, ib2)), eta);
			weight[ib2+ib1*GlobalV::NBANDS] = ( pelec->wg(ik,ib1)  - pelec->wg(iqk,ib2) )/factor/GlobalC::ucell.omega;
		}

	return;
}


void Chi0_standard:: Cal_last()
{
	for(int g0=0; g0<dim; g0++)
		for(int ib=0; ib<(oband*GlobalV::NBANDS); ib++)
		{
			B[ib][g0] = conj(b[g0][ib]);
		}

	return;
}

void Chi0_standard:: Cal_first()
{
	for(int g0=0; g0<dim; g0++)
		for(int ib=0; ib<(oband*GlobalV::NBANDS); ib++)
		{
			A[g0][ib] = weight[ib] * b[g0][ib];
		}

	return;
}

void Chi0_standard:: Cal_chi0(int iq, double omega, const elecstate::ElecState* pelec)
{
	for(int g0=0; g0<dim; g0++)
		for(int g1=0; g1<dim; g1++)
		{
			chi0[g0][g1] = std::complex<double>(0.0,0.0);
		}

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		int iqk = Cal_iq(ik, iq, GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
		Cal_b(iq, ik, iqk,GlobalC::rhopw);
		Cal_weight(iq, ik, omega, pelec);
		Cal_last();
		Cal_first();
		std::vector<std::vector<std::complex<double>>> A1(dim, std::vector<std::complex<double>>(oband*GlobalV::NBANDS));       // Peize Lin change ptr to std::vector at 2020.01.31
		std::vector<std::vector<std::complex<double>>> B1(oband*GlobalV::NBANDS, std::vector<std::complex<double>>(dim));       // Peize Lin change ptr to std::vector at 2020.01.31
		std::vector<std::vector<std::complex<double>>> C(dim, std::vector<std::complex<double>>(dim));                 // Peize Lin change ptr to std::vector at 2020.01.31

		for(int g0=0; g0<dim; g0++)
			for(int ib=0; ib<(oband*GlobalV::NBANDS); ib++)
			{
				A1[g0][ib] = A[g0][ib];
			}

		for(int ib=0; ib<(oband*GlobalV::NBANDS); ib++)
			for(int g0=0; g0<dim; g0++)
			{
				B1[ib][g0] = B[ib][g0];
			}

		std::complex<double> alpha=1;
		std::complex<double> beta=0;
		char transa = 'N';
		char transb = 'N';
		int M = dim;
		int N = dim;
		int K = oband*GlobalV::NBANDS;
		int lda = dim;
		int ldb = oband*GlobalV::NBANDS;
		int ldc = dim;
		zgemm_(&transa, &transb, &M, &N, &K, &alpha, ModuleBase::GlobalFunc::VECTOR_TO_PTR(B1[0]), &lda, ModuleBase::GlobalFunc::VECTOR_TO_PTR(A1[0]), &ldb, &beta, ModuleBase::GlobalFunc::VECTOR_TO_PTR(C[0]), &ldc);
		
		for(int g0=0; g0<dim; g0++)
			for(int g1=0; g1<dim; g1++)
			{
				chi0[g0][g1] += C[g0][g1];
			}
	}

	return;
}

void Chi0_standard:: Cal_rpa(int iq)
{
	for(int g0=0; g0<dim; g0++)
		for(int g1=0; g1<dim; g1++)
		{
			if(g0!=g1)
			{
				rpa[g0][g1] = -8.0 * ModuleBase::PI/qg2(iq,g0) * chi0[g0][g1];
			}
			else
			{
				rpa[g0][g1] = 1.0 - 8.0 * ModuleBase::PI/qg2(iq,g0) * chi0[g0][g1];			 	
			}
		}

	int l = Cinv(dim, rpa);
	if(l == 0)
	{
		ModuleBase::WARNING_QUIT("chi0_standard","(I-v*chi0) is a singular matrix !!");
	}

	return;
}

void Chi0_standard:: Cal_chi()
{
	CMatrixMul(dim,dim,dim,chi0,rpa,chi);

	return;
}

double Chi0_standard::qg2( int iq, int g0)
{
	double qg2;
	qg2 = ((GlobalC::kv.kvec_c[iq]+all_gcar[g0])*(ModuleBase::TWO_PI/GlobalC::ucell.lat0)).norm2();
	
	return qg2;
}
int Chi0_standard::Cinv(int n, std::complex<double>** a)
{
	int i,j,k,l,u,v,w;
	double p,q,s,d,b,m;
	std::complex<double> t;
	int *is, *js;
	is = new int[n];
	js = new int[n];

	for (k=0; k<=n-1; k++)
	{
		d=0.0;
		for (i=k; i<=n-1; i++)
			for (j=k; j<=n-1; j++)
			{
				p=a[i][j].real() * a[i][j].real() + a[i][j].imag() * a[i][j].imag();
				if (p>d)
				{
					d=p; is[k]=i; js[k]=j;
				}
			}

		if (d + 1.0 ==1.0)
		{
			delete[] is;
			delete[] js;
			std::cout << "error not inv" << std::endl;
			return(0);
		}

		if (is[k]!=k)
		{
			for (j=0; j<=n-1; j++)
			{
				t=a[k][j]; a[k][j]= a[is[k]][j]; a[is[k]][j]=t;
			}
		}

		if (js[k]!=k)
		{
			for (i=0; i<=n-1; i++)
			{
				t=a[i][k]; a[i][k]=a[i][js[k]]; a[i][js[k]]=t;
			}
		}

		l=k*n+k;
		a[k][k] = std::complex<double>( a[k][k].real()/d, -a[k][k].imag()/d);

		for (j=0; j<=n-1; j++)
		{
			if (j!=k)
			{
				p=a[k][j].real() *a[k][k].real(); q=a[k][j].imag() *a[k][k].imag();
				s=(a[k][j].real() +a[k][j].imag())*(a[k][k].real() + a[k][k].imag());
				a[k][j] = std::complex<double>( p-q, s-p-q);
			}
		}

		for (i=0; i<=n-1; i++)
		{
			if (i!=k)
			{
				for (j=0; j<=n-1; j++)
				{
					if (j!=k)
					{
						p=a[k][j].real() * a[i][k].real(); q=a[k][j].imag() *a[i][k].imag();
						s=(a[k][j].real() +a[k][j].imag())*(a[i][k].real() +a[i][k].imag());
						m=p-q; b=s-p-q;
						a[i][j] = std::complex<double>( a[i][j].real() - m, a[i][j].imag() - b);
					}
				}
			}
		}

		for (i=0; i<=n-1; i++)
		{
			if (i!=k)
			{
				p=a[i][k].real() *a[k][k].real(); q=a[i][k].imag() *a[k][k].imag();
				s=(a[i][k].real() +a[i][k].imag())*(a[k][k].real() +a[k][k].imag());
				a[i][k] = std::complex<double>( q-p, p+q-s);
			}
		}
	}

	for (k=n-1; k>=0; k--)
	{
		if (js[k]!=k)
		{
			for (j=0; j<=n-1; j++)
			{
				t=a[k][j]; a[k][j]=a[js[k]][j]; a[js[k]][j]=t;
			}
		}

		if (is[k]!=k)
		{
			for (i=0; i<=n-1; i++)
			{
				t=a[i][k]; a[i][k]= a[i][is[k]]; a[i][is[k]]=t;
			}
		}
	}

	delete[] is;
	delete[] js;

	return(1);
}


void Chi0_standard::CMatrixMul(int m, int n, int l, std::complex<double>** A, std::complex<double>** B, std::complex<double>** C)
{
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<l; j++)
		{
			C[i][j] = std::complex<double>(0.0,0.0);
			for(int k=0; k<n; k++)
			{
				C[i][j] = C[i][j] + A[i][k] * B[k][j] ;
			}
		}
	}

	return;
}

int Chi0_standard::Cal_iq(int ik, int iq, int a, int b, int c)
{
	int ik_x = ik%a;
	int ik_y = (ik%(a*b))/a;
	int ik_z = ik/(a*b);
	int iq_x = iq%a;
	int iq_y = (iq%(a*b))/a;
	int iq_z = iq/(a*b);
	int x = (ik_x + iq_x)%a;
	int y = (ik_y + iq_y)%b;
	int z = (ik_z + iq_z)%c;

	int Q = x + y * a + z * a * b;
	return Q;

}

int Chi0_standard::parallel_g()
{
	flag1 = new int[dim];  ModuleBase::GlobalFunc::ZEROS(flag1, dim);
	para_g = new double*[dim];
	for(int i=0;i<dim;i++)
	{
		para_g[i]= new double[2];
	}

	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<2;j++)
		{
			para_g[i][j] = -10000.0;
		}
	}

	bool exist;
	int num = 0;
	for(int i=0;i<dim;i++)
	{
		exist = false;
		for(int j=0;j<dim;j++)
		{
			if(all_gcar[i].x == para_g[j][0] && all_gcar[i].y == para_g[j][1])
			{
				flag1[i] = j;
				exist = true;
			}
		}
		if(!exist)
		{
			para_g[num][0] = all_gcar[i].x; para_g[num][1] = all_gcar[i].y;
			flag1[i] = num;
			num++;
		}
	}

	int dim1 = num;
	//std::cout << "dim1 = "<<dim1<<std::endl;
	for(int i=0;i<dim;i++)
	{
		std::cout <<"G["<<i<<"] = "<<all_gcar[i].x<<" "<<all_gcar[i].y<<" "<<all_gcar[i].z<<std::endl;
		std::cout <<"G_direct["<<i<<"] = "<<GlobalC::rhopw->gdirect[i].x<<" "<<GlobalC::rhopw->gdirect[i].y<<" "<<GlobalC::rhopw->gdirect[i].z<<std::endl;
		std::cout <<"flag1["<<i<<"] = "<<flag1[i]<<std::endl;
		std::cout <<"G_para["<<i<<"] = "<<para_g[i][0]<<"  "<<para_g[i][1]<<std::endl;
	}

	return dim1;
}

void Chi0_standard::chi0_para_g()
{
	for(int g0=0; g0<dim_para; g0++)
	{
		ModuleBase::GlobalFunc::ZEROS( chi0_para[g0], dim_para);
	}

	for(int g0=0; g0<dim; g0++)
	{
		for(int g1=0; g1<dim; g1++)
		{
			chi0_para[flag1[g0]][flag1[g1]] += chi0[g0][g1];
		}
	}

	return;
}
