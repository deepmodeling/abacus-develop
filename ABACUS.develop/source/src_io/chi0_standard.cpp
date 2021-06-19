//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-30
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating some optical properties using the linear response
//theory . 
//-------------------------------------------------------------------------------

#include "src_pw/global.h"
#include "chi0_standard.h"
#include "src_pw/hamilt_pw.h"
#include "src_pw/wavefunc_in_pw.h"
#include "src_io/optical.h"
#include "src_pw/klist.h"
#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

Chi0_standard chi0_standard;

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

void Chi0_standard:: Chi()
{
	TITLE("Chi0_standard","Chi");

	//---------------------------------------
	//  the serial number of q
	//---------------------------------------
	bool exist_q = false;
	for(int ik=0; ik<kv.nks; ik++)
	{
		double dx = fabs( kv.kvec_d[ik].x - q_start[0]);
		double dy = fabs( kv.kvec_d[ik].y - q_start[1]);
		double dz = fabs( kv.kvec_d[ik].z - q_start[2]);
		if( dx<0.0001 && dy<0.0001 && dz<0.0001)
		{
			start_q = ik; exist_q = true;
			break;
		}			
	}
	
	if(!exist_q)
	{
		WARNING_QUIT("chi0_hilbert","the chosen q is not included in the kmesh!!");
	}
	else
	{
		cout <<"start_q = "<<start_q<<endl;
	}	

	
	int icount = 0;
	std::vector<int> qc(kv.nks);		// Peize Lin change ptr to vector at 2020.01.31
	std::vector<double> ql(kv.nks);		// Peize Lin change ptr to vector at 2020.01.31
	int total_icount=0;
	int temp1; double temp2;
	
	if( direct[0]!=0 && direct[1]!=0 && direct[2]!=0)
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
			double p0 = x/direct[0]; double p1 = y/direct[1]; double p2 = z/direct[2];
			if( p0>0.0001 && fabs(p0-p1)<0.0001 && fabs(p0-p2)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}
	}
	else if( direct[0]==0 && direct[1]!=0 && direct[2]!=0 )
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
			double p1 = y/direct[1]; double p2 = z/direct[2];
			if( fabs(x)<0.0001 && p1>0.0001 && fabs(p1-p2)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p1; icount++;
			}
		}		
	}
	else if( direct[0]!=0 && direct[1]==0 && direct[2]!=0 )
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
			double p0 = x/direct[0]; double p2 = z/direct[2];
			if( fabs(y)<0.0001 && p0>0.0001 && fabs(p0-p2)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}		
	}		
	else if( direct[0]!=0 && direct[1]!=0 && direct[2]==0 )
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
			double p0 = x/direct[0]; double p1 = y/direct[1];
			if( fabs(z)<0.0001 && p0>0.0001 && fabs(p0-p1)<0.0001)
			{
				qc[icount] = ik; ql[icount] = p0; icount++;
			}
		}		
	}
	else if( direct[0]==0 && direct[1]==0 && direct[2]!=0 )
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
			double p2 = z/direct[2];
			if( fabs(x)<0.0001 && fabs(y)<0.0001 && p2 > 0.0001)
			{
				qc[icount] = ik; ql[icount] = p2; icount++;
			}
		}		
	}
	else if( direct[0]==0 && direct[1]!=0 && direct[2]==0 )
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
			double p1 = y/direct[1];
			if( fabs(x)<0.0001 && fabs(z)<0.0001 && p1 > 0.0001)
			{
				qc[icount] = ik; ql[icount] = p1; icount++;
			}
		}		
	}
	else if( direct[0]!=0 && direct[1]==0 && direct[2]==0 )
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			double x = kv.kvec_d[ik].x - kv.kvec_d[start_q].x;
			double y = kv.kvec_d[ik].y - kv.kvec_d[start_q].y;
			double z = kv.kvec_d[ik].z - kv.kvec_d[start_q].z;
			
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
	//	WARNING_QUIT("chi0_hilbert","Now the kmesh contains no kpoint along this direction!");
	//}
	if(total_icount < nq-1 )
	{
		WARNING_QUIT("chi0_hilbert","Now the kmesh doesn't contain enough kpoints along this direction! please change the parameter nq smaller");
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
	//cout <<"temp_q = "<<temp_q<<endl;
	
	interval_q = temp_q - start_q;
	//cout <<"interval_q = "<<interval_q<<endl;
	
	Parallel_G();
	
	if(system == "surface")
	{
		dim_para = parallel_g();
		cout << "dim_para = "<<dim_para<<endl;
	}
	
	//----------------------------------------------------------
	// calculate the number of occupied bands
	//----------------------------------------------------------	
	int occ_bands = 0;
	bool occ_flag;
	for(int ib=0; ib<NBANDS; ib++)
	{
		occ_flag = false;
		for(int ik=0; ik<kv.nks; ik++)
		{
			if( wf.wg(ik,ib)> 0.0001)
			{
				occ_flag = true;
				continue;
			}
		}
		
		if(occ_flag == true)
			occ_bands++;
	}
	cout <<"occ_bands = "<<occ_bands<<endl;	
	oband = occ_bands;
	
	if(!init_finish)
	{
		Init();
		init_finish = true;
	}
	
	for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
	{
		double q = sqrt(((kv.kvec_c[iq])*(TWO_PI/ucell.lat0)).norm2());
		double gather[nomega];
		
		int count =0;
		for(double omega=0.0; omega<(domega*nomega); omega=omega+domega)
		{
			Cal_chi0(iq,omega);
			cout<<"chi0 iq= "<<iq<<" omega= "<<omega<<"  "<<chi0[0][0].real()<<" "<<chi0[0][0].imag()<<endl;
			if(system == "surface")
			{
				chi0_para_g();
				cout<<"chi0_para iq= "<<iq<<" omega= "<<omega<<"  "<<chi0_para[0][0].real()<<" "<<chi0_para[0][0].imag()<<endl;
			}						
			Cal_rpa(iq);
			Cal_chi();
			cout<<"chi iq= "<<iq<<" omega= "<<omega<<"  "<<chi[0][0].real()<<" "<<chi[0][0].imag()<<endl;
			cout<<"epsilon iq= "<<iq<<" omega= "<<omega<<"  "<<8*PI/q/q*chi[0][0].real()<<" "<<8*PI/q/q*chi[0][0].imag()<<endl;
			gather[count] = -8*PI/q/q*chi[0][0].imag();		
			count++;
		}

		if(out_epsilon)
		{
			stringstream sseps;
			sseps << global_out_dir << "Imeps^-1_"<<iq<<".dat";
			ofstream ofseps(sseps.str().c_str());
			ofseps<<"Energy(Ry)"<<"   "<<"-Im{epsilon^-1}"<<endl;
			for(int i=0; i<nomega; i++)
			{
				ofseps <<i*domega<<"   "<<gather[i]<<endl;
			}
			ofseps.close();
		}		
	}
	
	Delete();
	return;
}

void Chi0_standard::Parallel_G()
{
	TITLE("Chi0_standard","Parallel_G");
	//----------------------------
	// init
	//----------------------------
	num_G_core = new int[DSIZE];
	num_G_dis = new int[DSIZE];
	G_r_core = new double[pw.ngmc];
	num_Gvector_core = new int[DSIZE];
	num_Gvector_dis = new int[DSIZE];
	G_r = new double[pw.ngmc_g];
	Gvec_core = new double[3*pw.ngmc];
	Gvec = new double[3*pw.ngmc_g];
	all_gcar = new Vector3<double>[pw.ngmc_g];
	flag = new int[pw.ngmc_g];
	
	for(int i=0;i<pw.ngmc_g;i++)
	{
		flag[i] = i;
	}
	
	ZEROS( num_G_dis, DSIZE);
	ZEROS( num_G_core, DSIZE);
	ZEROS( num_Gvector_dis, DSIZE);
	ZEROS( num_Gvector_core, DSIZE);
	
#ifdef __MPI
	MPI_Allgather( &pw.ngmc, 1, MPI_INT, num_G_core, 1, MPI_INT, POOL_WORLD);
#endif
	
	memset(num_G_dis,0,DSIZE*sizeof(int));
	for(int i=0; i<DSIZE; i++)
	{
		for(int j=0; j<i; j++)
		{
			num_G_dis[i] += num_G_core[j];
		}
	}
	
	for(int i=0;i<DSIZE;i++)
	{
		num_Gvector_dis[i] = num_G_dis[i] * 3;
		num_Gvector_core[i] = num_G_core[i] * 3;
	}
	
	for(int g0=0;g0<pw.ngmc; g0++)
	{
		G_r_core[g0] = pw.get_NormG_cartesian(g0);
		Gvec_core[3 * g0] = pw.get_G_cartesian_projection(g0, 0);
		Gvec_core[3 * g0 + 1] = pw.get_G_cartesian_projection(g0, 1);
		Gvec_core[3 * g0 + 2] = pw.get_G_cartesian_projection(g0, 2);
	}
	
#ifdef __MPI
	MPI_Allgatherv( G_r_core, pw.ngmc, MPI_DOUBLE, G_r, num_G_core, num_G_dis, MPI_DOUBLE, POOL_WORLD);
	MPI_Allgatherv( Gvec_core, 3*pw.ngmc, MPI_DOUBLE, Gvec, num_Gvector_core, num_Gvector_dis, MPI_DOUBLE, POOL_WORLD);
#endif
	
	double t1; int t2;
	for(int i=0;i<pw.ngmc_g;i++)
	{
		for(int j=0;j<pw.ngmc_g-i-1;j++)
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
	
	for(int i=0;i<pw.ngmc_g;i++)
	{
		all_gcar[i].x = Gvec[3*i]; all_gcar[i].y = Gvec[3*i+1]; all_gcar[i].z = Gvec[3*i+2];
		//cout<<"all_gcar["<<i<<"]= "<<all_gcar[i].x<<" "<<all_gcar[i].y<<" "<<all_gcar[i].z<<endl;
	}
	
	return;
}

void Chi0_standard:: Init()
{

	b_core = new complex<double>[pw.ngmc];  
	
	b_summary = new complex<double>[pw.ngmc_g]; 
	
	b_order = new complex<double>[pw.ngmc_g];

	psi_r1 = new complex<double>*[NBANDS];
	for(int ib=0; ib<NBANDS; ib++)
	{
		psi_r1[ib] = new complex<double>[pw.nrxx];
	}
		
	psi_r2 = new complex<double>*[NBANDS];
	for(int ib=0; ib<NBANDS; ib++)
	{
		psi_r2[ib] = new complex<double>[pw.nrxx];
	}	
	cout << "psi OK" <<endl;
	
	b = new complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		b[g0] = new complex<double>[oband*NBANDS];
	}
	cout << "b ok"<<endl;
	
	A = new complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		A[g0] = new complex<double>[oband*NBANDS];
	}
	cout << "A ok"<<endl;	
	
	B = new complex<double>*[oband*NBANDS];
	for(int ib=0; ib<(oband*NBANDS); ib++)
	{
		B[ib] = new complex<double>[dim];
	}
	cout << "B ok" << endl;
	
	weight = new complex<double> [oband*NBANDS];
	cout << "weight OK"<<endl;
	
	chi0 = new complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		chi0[g0] = new complex<double>[dim];
	}
	cout << "chi0 OK" <<endl;
	
	chi0_para = new complex<double>*[dim_para];
	for(int g0=0; g0<dim_para; g0++)
	{
		chi0_para[g0] = new complex<double>[dim_para];
	}

	rpa = new complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		rpa[g0] = new complex<double>[dim];
	}
	cout << "rpa OK" <<endl;

	chi = new complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		chi[g0] = new complex<double>[dim];
	}
	cout << "chi OK" <<endl;	
	
	cout << "All ok"<<endl;
		
	return;
}

void Chi0_standard::Delete()
{
	TITLE("Chi0_standard","Delete");
	if(init_finish == true)				// Peize Lin change = to == at 2020.01.31
	{
		delete[] b_core;
		delete[] b_summary;
		delete[] b_order;
		
		for(int ib=0; ib<NBANDS; ib++)
		{
			delete[] psi_r1[ib];
		}
		delete[] psi_r1;
			
		for(int ib=0; ib<NBANDS; ib++)
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

		for(int ib=0; ib<(oband*NBANDS); ib++)
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

void Chi0_standard::Cal_Psi(int iq, complex<double> **psi_r)
{
	double phase_x, phase_xy, phase_xyz;
	complex<double> exp_tmp;
	for(int ib = 0; ib < NBANDS; ib++)
	{
		ZEROS( UFFT.porter, (pw.nrxx) );
		for(int ig = 0; ig < kv.ngk[iq] ; ig++)
		{
			UFFT.porter[ pw.ig2fftw[wf.igk(iq,ig)] ] = wf.evc[iq](ib,ig);
		}
		
		pw.FFT_wfc.FFT3D(UFFT.porter,1);
		int ir=0;
		for(int ix=0; ix<pw.ncx; ix++)
		{
			phase_x = kv.kvec_d[iq].x*ix/pw.ncx;
			for(int iy=0; iy<pw.ncy; iy++)
			{
				phase_xy = phase_x + kv.kvec_d[iq].y*iy/pw.ncy;
				for(int iz=pw.nczp_start; iz<pw.nczp_start+pw.nczp; iz++)
				{
					phase_xyz = (phase_xy + kv.kvec_d[iq].z*iz/pw.ncz) *TWO_PI;
					exp_tmp = complex<double>( cos(phase_xyz), sin(phase_xyz) );
					psi_r[ib][ir] = UFFT.porter[ir]*exp_tmp;
					ir++;
				}
					
			}
		}
	}
	
	return;
}

void Chi0_standard::Cal_b(int iq, int ik, int iqk)
{
	TITLE("Chi0_standard","Cal_b");
	Vector3<double> qk;
	qk = kv.kvec_c[iq] + kv.kvec_c[ik];
	//cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<endl;
	double phase_x, phase_xy, phase_xyz;
	Vector3<double> q = kv.kvec_d[iq];
	complex<double> exp_tmp;
	
	Cal_Psi(ik, psi_r1);
	Cal_Psi(iqk, psi_r2);
	
	for(int ib1=0; ib1< oband; ib1++)
	{
		for(int ib2=0; ib2<NBANDS; ib2++)
		{
			int ir=0;
			for(int ix=0; ix<pw.ncx; ix++)
			{
				phase_x = q.x*ix/pw.ncx;
				for(int iy=0; iy<pw.ncy; iy++)
				{
					phase_xy = phase_x + q.y*iy/pw.ncy;
					for(int iz=pw.nczp_start; iz<pw.nczp_start+pw.nczp; iz++)
					{
						phase_xyz = (phase_xy + q.z*iz/pw.ncz) *TWO_PI;
						exp_tmp = complex<double>(cos(-phase_xyz), sin(-phase_xyz));
						UFFT.porter[ir] = conj(psi_r1[ib1][ir]) * psi_r2[ib2][ir] *exp_tmp;
						ir++;
					}
				}
			}
			
			pw.FFT_chg.FFT3D( UFFT.porter, -1);
			//for(int g0=0; g0<dim; g0++)
			//{
			//	b[g0][ib1][ib2] = UFFT.porter[ pw.ig2fftc[g0] ];
			//}
			for(int g0=0;g0<pw.ngmc; g0++)
			{
				b_core[g0] = UFFT.porter[ pw.ig2fftc[g0] ];
			}
			
#ifdef __MPI
			MPI_Allgatherv( b_core, pw.ngmc, mpicomplex, b_summary, num_G_core, num_G_dis, mpicomplex, POOL_WORLD);
#endif
			
			for(int i=0;i<pw.ngmc_g;i++)
			{
				b_order[i] = b_summary[flag[i]];
			}
			
			for(int g0=0; g0<dim; g0++)
			{
				b[g0][ib2+ib1*NBANDS] = b_order[g0];
			}

		}
	}
	
	return;
}

void Chi0_standard:: Cal_weight(int iq, int ik, double omega)
{
	int iqk = Cal_iq(ik, iq, kv.nmp[0], kv.nmp[1], kv.nmp[2]);
	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<NBANDS; ib2++)
		{
			complex<double> factor = complex<double>( (omega + wf.ekb[ik][ib1] - wf.ekb[iqk][ib2]), eta);
			weight[ib2+ib1*NBANDS] = ( wf.wg(ik,ib1)  - wf.wg(iqk,ib2) )/factor/ucell.omega;
		}
		
	return;
}


void Chi0_standard:: Cal_last()
{
	for(int g0=0; g0<dim; g0++)
		for(int ib=0; ib<(oband*NBANDS); ib++)
		{
			B[ib][g0] = conj(b[g0][ib]);
		}
		
	return;
}

void Chi0_standard:: Cal_first()
{
	for(int g0=0; g0<dim; g0++)
		for(int ib=0; ib<(oband*NBANDS); ib++)
		{
			A[g0][ib] = weight[ib] * b[g0][ib];
		}
		
	return;
}

void Chi0_standard:: Cal_chi0(int iq, double omega)
{
	for(int g0=0; g0<dim; g0++)
		for(int g1=0; g1<dim; g1++)
		{
			chi0[g0][g1] = complex<double>(0.0,0.0);
		}
	
	for(int ik=0; ik<kv.nks; ik++)
	{
		int iqk = Cal_iq(ik, iq, kv.nmp[0], kv.nmp[1], kv.nmp[2]);
		Cal_b(iq, ik, iqk);
		Cal_weight(iq, ik, omega);
		Cal_last();
		Cal_first();
		std::vector<std::vector<complex<double>>> A1(dim, std::vector<complex<double>>(oband*NBANDS));       // Peize Lin change ptr to vector at 2020.01.31
		std::vector<std::vector<complex<double>>> B1(oband*NBANDS, std::vector<complex<double>>(dim));       // Peize Lin change ptr to vector at 2020.01.31
		std::vector<std::vector<complex<double>>> C(dim, std::vector<complex<double>>(dim));                 // Peize Lin change ptr to vector at 2020.01.31
		
		for(int g0=0; g0<dim; g0++)
			for(int ib=0; ib<(oband*NBANDS); ib++)
			{
				A1[g0][ib] = A[g0][ib];
			}
		
		for(int ib=0; ib<(oband*NBANDS); ib++)
			for(int g0=0; g0<dim; g0++)
			{
				B1[ib][g0] = B[ib][g0];
			}
			
		complex<double> alpha=1;
		complex<double> beta=0;
		char transa = 'N';
		char transb = 'N';
		int M = dim;
		int N = dim;
		int K = oband*NBANDS;
		int lda = dim;
		int ldb = oband*NBANDS;
		int ldc = dim;
		zgemm_(&transa, &transb, &M, &N, &K, &alpha, VECTOR_TO_PTR(B1[0]), &lda, VECTOR_TO_PTR(A1[0]), &ldb, &beta, VECTOR_TO_PTR(C[0]), &ldc);
		
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
				rpa[g0][g1] = -8.0 * PI/qg2(iq,g0) * chi0[g0][g1];
			}
			else
			{
				rpa[g0][g1] = 1.0 - 8.0 * PI/qg2(iq,g0) * chi0[g0][g1];			 	
			}
		}
		
	int l = Cinv(dim, rpa);
	if(l == 0)
	{
		WARNING_QUIT("chi0_standard","(I-v*chi0) is a singular matrix !!");
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
	qg2 = ((kv.kvec_c[iq]+all_gcar[g0])*(TWO_PI/ucell.lat0)).norm2();
	
	return qg2;
}
int Chi0_standard::Cinv(int n, complex<double>** a)
{
	int i,j,k,l,u,v,w;
	double p,q,s,d,b,m;
	complex<double> t;
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
			cout << "error not inv" << endl;
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
		a[k][k] = complex<double>( a[k][k].real()/d, -a[k][k].imag()/d);
		
		for (j=0; j<=n-1; j++)
		{
			if (j!=k)
			{
				p=a[k][j].real() *a[k][k].real(); q=a[k][j].imag() *a[k][k].imag();
				s=(a[k][j].real() +a[k][j].imag())*(a[k][k].real() + a[k][k].imag());
				a[k][j] = complex<double>( p-q, s-p-q);
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
						a[i][j] = complex<double>( a[i][j].real() - m, a[i][j].imag() - b);
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
				a[i][k] = complex<double>( q-p, p+q-s); 
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


void Chi0_standard::CMatrixMul(int m, int n, int l, complex<double>** A, complex<double>** B, complex<double>** C)
{
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<l; j++)
		{
			C[i][j] = complex<double>(0.0,0.0);
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
	flag1 = new int[dim];  ZEROS(flag1, dim);
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
	//cout << "dim1 = "<<dim1<<endl;
	for(int i=0;i<dim;i++)
	{
		cout <<"G["<<i<<"] = "<<all_gcar[i].x<<" "<<all_gcar[i].y<<" "<<all_gcar[i].z<<endl;
		cout <<"G_direct["<<i<<"] = "<<pw.gdirect[i].x<<" "<<pw.gdirect[i].y<<" "<<pw.gdirect[i].z<<endl;
		cout <<"flag1["<<i<<"] = "<<flag1[i]<<endl;
		cout <<"G_para["<<i<<"] = "<<para_g[i][0]<<"  "<<para_g[i][1]<<endl;  
	}
	
	return dim1;
}

void Chi0_standard::chi0_para_g()
{
	for(int g0=0; g0<dim_para; g0++)
	{
		ZEROS( chi0_para[g0], dim_para);
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



