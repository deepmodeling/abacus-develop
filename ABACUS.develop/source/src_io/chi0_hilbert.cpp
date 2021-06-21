//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-14
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating some optical properties using the linear response
//theory . and it is a basic work and use the hilbert-transform to reduce the
//calculated quantities.
//-------------------------------------------------------------------------------
#include "src_pw/global.h"
#include "chi0_hilbert.h"
#include "src_pw/hamilt_pw.h"
#include "src_pw/wavefunc_in_pw.h"
#include "../src_io/optical.h"
#include "src_pw/klist.h"
#include <iostream>
#include "../src_lcao/global_fp.h"

using namespace std;

Chi0_hilbert chi0_hilbert;

Chi0_hilbert::Chi0_hilbert()
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
	out_chi = false;           
	out_chi0 = false;          
	fermi_level = 0.0;     
	coulomb_cutoff = false;    
	kmesh_interpolation = false;	
	for(int i=0; i<100; i++)
	{
		qcar[i][0] = 0.0; qcar[i][1] = 0.0; qcar[i][2] = 0.0;
	}	
	lcao_box[0] = 10; lcao_box[1] = 10; lcao_box[2] = 10;
}

Chi0_hilbert::~Chi0_hilbert()
{
}

// begin calculate
#include "src_pw/occupy.h"
void Chi0_hilbert::Chi()
{
	TITLE("Chi0_hilbert","Chi");
	//---------------------------------------
	// check nomega
	//---------------------------------------
	double energy_scale = domega * nomega;
	cout <<"energy_scale = "<<energy_scale<<endl;
	double max_e = 0.0;
	for(int ik1=0; ik1<kv.nks; ik1++)
	{
		for(int ik2=0; ik2<kv.nks; ik2++)
		{
			if(max_e < (wf.ekb[ik1][NBANDS-1] - wf.ekb[ik2][0]))
			{
				max_e = wf.ekb[ik1][NBANDS-1] - wf.ekb[ik2][0];
			}			
		}		 
	}
	cout <<"max_e = "<<max_e<<endl;
	if(max_e > energy_scale)
	{
		int nomega_right = int(max_e/domega);
		cout << "nomega must >= "<<nomega_right<<endl;
		WARNING_QUIT("chi0_hilbert","nomega is too small!!!!!!!!!");
	}
	
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
	std::vector<int> qc(kv.nks); 			// Peize Lin change ptr to vector at 2020.01.31
	std::vector<double> ql(kv.nks);			// Peize Lin change ptr to vector at 2020.01.31
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
	
	//----------------------------------------------------------
	//    parallel in G     collection and sorting
	//----------------------------------------------------------
	Parallel_G();
	
	if(system == "surface")    // calculate the g_//
	{
		dim_para = parallel_g();
		//cout << "dim_para = "<<dim_para<<endl;
	}

	//cout <<"en.ef = "<<en.ef<<endl;
	//cout <<"degauss = "<<INPUT.degauss<<endl;
	
	//------------------------------------------------------------------
	// change the fermi level or not 
	//------------------------------------------------------------------
	double energy = en.ef + fermi_level;
	//cout <<"energy = "<<energy<<endl;
	
	cweight = new double*[kv.nks]; 
	for(int i=0;i<kv.nks;i++)             // new occupation
	{
		cweight[i] = new double[NBANDS];
	}
	
	for(int ik = 0;ik < kv.nks;ik++)
	{
		for(int ib=0; ib<NBANDS; ib++)
		{
			//----------------------------------------------------------
			//  gauss smearing
			//----------------------------------------------------------
			if(INPUT.smearing == "gauss" || INPUT.smearing == "gaussian")
			{
				cweight[ik][ib] = kv.wk[ik] * Occupy::wgauss( (energy - wf.ekb[ik][ib] )/INPUT.degauss, 0);
			}									
			//----------------------------------------------------------
			//  fixed smearing
			//----------------------------------------------------------
			else if(INPUT.smearing == "fixed")
			{
				if((wf.ekb[ik][ib]-energy)<-0.0001 )
				{
					cweight[ik][ib] = kv.wk[ik] * 1.0; 
				}
				else if(fabs(wf.ekb[ik][ib]-energy)<0.0001 )
				{
					cweight[ik][ib] = kv.wk[ik] * 0.5;
				}
				else
				{
					cweight[ik][ib] = 0.0;
				}
			}						
			//----------------------------------------------------------
			//   fermi-dirac smearing
			//----------------------------------------------------------
			else if(INPUT.smearing == "fd")
			{
				cweight[ik][ib] = kv.wk[ik] * Occupy::wgauss( (energy - wf.ekb[ik][ib] )/INPUT.degauss, -99);
			}
			else
			{
				WARNING_QUIT("chi0_hilbert","Sorry! This smearing hasn't been realized here yet!!!");
			}			
		} 
	}
	
	//----------------------------------------------------------
	//  calculate the total electrons
	//----------------------------------------------------------
	double elec = 0.0;
	for(int ik = 0;ik < kv.nks;ik++)
	{
		for(int ib=0; ib<NBANDS; ib++)
		{
			elec += cweight[ik][ib];
		}
	}
	
	cout <<"elec = "<<elec<<endl;
	
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
			if( cweight[ik][ib] > 0.0001)
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
	
	//----------------------------------------------------------
	// init
	//----------------------------------------------------------
	
	if(!init_finish)
	{
		Init();
		init_finish = true;
	}
		
	//----------------------------------------------------------
	// bulk calculating
	//----------------------------------------------------------
	if(system == "bulk")
	{
		//-------------------------------------------------
		// pw method
		//-------------------------------------------------
		if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
		{
			for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
			{
				time_t start,end;
				start = time(NULL);
				Cal_Chi0s(iq);  // calculate chi0s
				end = time(NULL);
				cout<<"chi0s time:: "<<end-start<<" s"<<endl;
				
				time_t start1,end1;
				start1 = time(NULL);
				Cal_T();
				end1 = time(NULL);
				cout<<"T time:: "<<end1-start1<<" s"<<endl;
				
				time_t start2,end2;
				start2 = time(NULL);
				Cal_Chi0();   // calculate  chi0
				end2 = time(NULL);
				cout<<"chi0 time:: "<<end2-start2<<" s"<<endl;
				
				if(out_chi0)
				{
					plot_chi0(iq);
				}
				for(int i=0; i<nomega; i++)
				{
					cout <<"iq = "<<iq<<" chi0["<<i<<"] = "<<chi0[0][i].real()<<"  "<<chi0[0][i].imag()<<endl;
				}
				
				Cal_kernel(iq);
				
				time_t start3,end3;
				start3 = time(NULL);
				Cal_Chi(iq);  // calculate chi
				end3 = time(NULL);
				cout<<"chi time:: "<<end3-start3<<" s"<<endl;
			}					
		}
		//-------------------------------------------------
		// lcao method
		//-------------------------------------------------
		else if(BASIS_TYPE == "lcao")
		{
			//-------------------------------------------------
			// first, readin files( "nearest.dat" and "overlap_q")
			//-------------------------------------------------
			int tmp; int icount = 0;
			//------------------------
			// nearest.dat
			//------------------------
			stringstream ss;
			ss << global_out_dir <<"nearest.dat";
			cout << ss.str().c_str() << endl;
			ifstream ifsn(ss.str().c_str());
			if(!ifsn)
			{
				WARNING_QUIT("chi0_hilbert", "Can't find the nearest.dat file!");
			}
			ifsn >> NR;
			cout << "NR in = "<<NR << endl;
			
			R = new Vector3<int>** [NLOCAL];
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				R[iw1] = new Vector3<int>* [NLOCAL];
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					R[iw1][iw2] = new Vector3<int>[NR];
				}
			}
			
			Rcar = new Vector3<double>** [NLOCAL];
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				Rcar[iw1] = new Vector3<double>* [NLOCAL];
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					Rcar[iw1][iw2] = new Vector3<double>[NR];
				}
			}
			
			Rmax = new int* [NLOCAL];
			for(int iw=0; iw<NLOCAL; iw++)
			{
				Rmax[iw] = new int[NLOCAL];
			}
			
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					//cout <<"iw1 = "<<iw1<<" iw2 = "<<iw2<<endl;
					ifsn >> tmp; 
					ifsn >> tmp;
					ifsn >> Rmax[iw1][iw2]; //cout << Rmax[iw1][iw2] << endl;
					for(int i=0; i<Rmax[iw1][iw2]; i++)
					{
						ifsn >> R[iw1][iw2][i].x; //cout <<R[iw1][iw2][i].x;
						ifsn >> R[iw1][iw2][i].y; //cout <<R[iw1][iw2][i].y;
						ifsn >> R[iw1][iw2][i].z; //cout <<R[iw1][iw2][i].z; cout << endl;
					}
				}
			}
			ifsn.close();
			
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					for(int i=0; i<Rmax[iw1][iw2]; i++)
					{
						Rcar[iw1][iw2][i].x = ucell.latvec.e11 * R[iw1][iw2][i].x 
						+ ucell.latvec.e21 * R[iw1][iw2][i].y + ucell.latvec.e31 * R[iw1][iw2][i].z;
						Rcar[iw1][iw2][i].y = ucell.latvec.e12 * R[iw1][iw2][i].x 
						+ ucell.latvec.e22 * R[iw1][iw2][i].y + ucell.latvec.e32 * R[iw1][iw2][i].z;
						Rcar[iw1][iw2][i].z = ucell.latvec.e13 * R[iw1][iw2][i].x 
						+ ucell.latvec.e23 * R[iw1][iw2][i].y + ucell.latvec.e33 * R[iw1][iw2][i].z;
					}
				}
			}
			
			overlap = new complex<double> ***[NLOCAL];
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				overlap[iw1] = new complex<double>**[NLOCAL];
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					overlap[iw1][iw2] = new complex<double>*[dim];
					for(int ig=0; ig<dim; ig++)
					{
						overlap[iw1][iw2][ig] = new complex<double>[NR];
						ZEROS( overlap[iw1][iw2][ig], NR);
					}
				}
			}
			
			for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
			{
				//--------------------------------
				// readin overlap_q
				//--------------------------------
				//cout << "begin overlap"<<endl;
				stringstream ssq;
				ssq << global_out_dir <<"q_"<<icount;
				ifstream ifso(ssq.str().c_str());
				if (!ifso)
				{
					WARNING_QUIT("chi0_hilbert", "Can't find the overlap_q file!");
				}
				ifso >> tmp; //cout << tmp << endl;
				for(int ig=0; ig<dim; ig++)
				{
					for(int iw1=0; iw1<NLOCAL; iw1++)
					{
						for(int iw2=0; iw2<NLOCAL; iw2++)
						{
							for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
							{
								ifso >> overlap[iw1][iw2][ig][ir];
								//cout<<"overlap["<<iw1<<"]["<<iw2<<"]["<<ig<<"]["<<ir<<"] = "<<overlap[iw1][iw2][ig][ir] <<endl;
							}
						}
					}
				}
				ifso.close();
				
				time_t start,end;
				start = time(NULL);
				Cal_Chi0s(iq);  // calculate chi0s
				end = time(NULL);
				cout<<"chi0s time:: "<<end-start<<" s"<<endl;
				
				time_t start1,end1;
				start1 = time(NULL);
				Cal_T();
				end1 = time(NULL);
				cout<<"T time:: "<<end1-start1<<" s"<<endl;
				
				time_t start2,end2;
				start2 = time(NULL);
				Cal_Chi0();   // calculate  chi0
				end2 = time(NULL);
				cout<<"chi0 time:: "<<end2-start2<<" s"<<endl;
				
				if(out_chi0)
				{
					plot_chi0(iq);
				}
				for(int i=0; i<nomega; i++)
				{
					cout <<"iq = "<<iq<<" chi0["<<i<<"] = "<<chi0[0][i].real()<<"  "<<chi0[0][i].imag()<<endl;
				}
				
				Cal_kernel(iq);
				
				time_t start3,end3;
				start3 = time(NULL);
				Cal_Chi(iq);  // calculate chi
				end3 = time(NULL);
				cout<<"chi time:: "<<end3-start3<<" s"<<endl;
				
				icount++;
			}
			
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					delete[] R[iw1][iw2];
				}
				delete[] R[iw1];
			}
			delete[] R;
			
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					delete[] Rcar[iw1][iw2];
				}
				delete[] Rcar[iw1];
			}
			delete[] Rcar;
			
			for(int iw=0; iw<NLOCAL; iw++)
			{
				delete[] Rmax[iw];
			}
			delete[] Rmax;
			
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					for(int ig=0; ig<dim; ig++)
					{
						delete[] overlap[iw1][iw2][ig];
					}
					delete[] overlap[iw1][iw2];
				}
				delete[] overlap[iw1];
			}
			delete[] overlap;
												
		}
		else
		{
			WARNING_QUIT("chi0_hilbert","pw, lcao_in_pw, or lcao must be chosen!!!!!!");
		}
	}
	//----------------------------------------------------------
	// surface calculating( maybe some bugs in it !!!!!!!!!!)
	//----------------------------------------------------------
	else if(system == "surface")
	{
		for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
		{
			time_t start, end;
			start = time(NULL);
			Cal_Chi0s(iq);
			end = time(NULL);
			cout<<"chi0s time:: "<<end-start<<" s"<<endl;
			
			time_t start1, end1;
			start1 = time(NULL);
			Cal_T();
			end1 = time(NULL);
			cout<<"T time:: "<<end1-start1<<" s"<<endl;
			
			time_t start2, end2;
			start2 = time(NULL);
			Cal_Chi0();
			end2 = time(NULL);
			cout<<"chi time:: "<<end2-start2<<" s"<<endl;
			
			if(out_chi0)
			{
				plot_chi0(iq);
			}
			for(int i=0; i<nomega; i++)
			{
				cout <<"iq = "<<iq<<" chi0["<<i<<"] = "<<chi0[0][i].real()<<"  "<<chi0[0][i].imag()<<endl;
			}
			
			Cal_kernel_2D(iq);
			
			time_t start3, end3;
			start3 = time(NULL);
			Cal_Chi_surface(iq);
			end3 = time(NULL);
			cout<<"chi time:: "<<end3-start3<<" s"<<endl;
		}
	}
	else
	{
		WARNING_QUIT("chi0_hilbert","system must be bulk or surface!!!!!!!");
	}
	
	Delete();
	
	return; 	
}

void Chi0_hilbert::Init()
{
	TITLE("Chi0_hilbert","Init");
	//cout << "nbands(init) = " <<NBANDS <<endl;
	//cout << "oband = " <<oband <<endl;
	//cout << "nrxx = "<<pw.nrxx<<endl;
	//cout << "ngmc = " << pw.ngmc <<endl;
	//cout << "kv.ngk[0] = "<<kv.ngk[0]<<endl;
	//cout << "dim = "<<dim<<endl;
	//cout << "kv.nks = " <<kv.nks <<endl;
	//cout << "kv.nmp = " << kv.nmp[0]<<"  "<<kv.nmp[1]<<"  "<<kv.nmp[2]<<endl;
	
	b_core = new complex<double>[pw.ngmc];  
	
	b_summary = new complex<double>[pw.ngmc_g]; 
	
	b_order = new complex<double>[pw.ngmc_g];
	
	/*psi_r = new complex<double>**[kv.nks];
	for(int iq=0; iq<kv.nks; iq++)
	{
		psi_r[iq] = new complex<double>*[NBANDS];
		for(int ib=0; ib<NBANDS; ib++)
		{
			psi_r[iq][ib] = new complex<double>[pw.nrxx];
		}
	}*/
	
	if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
	{
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
	}
	//cout << "psi1 OK" <<endl;
	//cout << "psi2 OK" <<endl;
	cout << "psi OK" <<endl;
	
	b = new complex<double>**[dim];
	for(int g0=0; g0<dim; g0++)
	{
		b[g0] = new complex<double>*[oband];
		for(int ib=0; ib<oband; ib++)
		{
			b[g0][ib] = new complex<double>[NBANDS];
		}
	}
	cout << "b OK" <<endl;
	
	chi0s = new complex<double>*[dim*dim];       // calculate the chi0s (imagniary part of chi0)
	for(int g0=0; g0<dim*dim; g0++)
	{
		chi0s[g0] = new complex<double>[nomega];
	}
	cout << "chi0s OK" <<endl;
	
	chi0 = new complex<double>*[dim*dim];
	for(int g0=0; g0<dim*dim; g0++)
	{
		chi0[g0] = new complex<double>[nomega];
	}
	cout << "chi0 OK" <<endl;
	
	chi0_gg = new complex<double>*[dim];       // split the nomega of chi0
	for(int g0=0; g0<dim; g0++)
	{
		chi0_gg[g0] = new complex<double> [dim];
	}
	cout << "chi0_gg OK" <<endl;
	
	T = new complex<double>*[nomega];              // tranaslation  matrix
	for(int i=0; i<nomega; i++)
	{
		T[i] = new complex<double>[nomega];
	}
	cout << "T OK" <<endl;
		
	rpa = new complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		rpa[g0] = new complex<double> [dim];
	}
	cout << "rpa OK" <<endl;
	
	cout <<"system = "<<system<<endl;
	if(system == "surface")
	{
		chi0_para = new complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			chi0_para[g0] = new complex<double> [dim_para];
		}
		cout << "chi0_para OK" <<endl;
		
		kernel = new complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			kernel[g0] = new complex<double> [dim_para];
		}
		cout << "kernel OK" <<endl;

		rpa1 = new complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			rpa1[g0] = new complex<double> [dim_para];
		}
		cout << "rpa1 OK" <<endl;
		
		chi_para = new complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			chi_para[g0] = new complex<double> [dim_para];
		}
		cout << "chi_para" <<endl;
	}
	else
	{
		kernel = new complex<double>*[dim];
		for(int g0=0; g0<dim; g0++)
		{
			kernel[g0] = new complex<double> [dim];
		}
		cout << "kernel OK" <<endl;
		
		chi = new complex<double>*[dim];
		for(int g0=0; g0<dim; g0++)
		{
			chi[g0] = new complex<double> [dim];
		}
		cout << "chi OK" <<endl;
	}
	cout << "ALL OK"<<endl;
	
	return;	
}

void Chi0_hilbert::Delete()
{
	TITLE("Chi0_hilbert","Delete");
	if(init_finish)
	{
		delete[] b_core;
		delete[] b_summary;
		delete[] b_order;
		
		for(int i=0; i<kv.nks; i++)
		{
			delete[] cweight[i];
		}
		delete[] cweight;
		
		if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
		{
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
		}
		
		for(int g0=0;g0<dim; g0++)
		{
			for(int ib=0; ib<oband; ib++)
			{
				delete[] b[g0][ib];
			}
			delete[] b[g0];
		}
		delete[] b;
		
		for(int g0=0;g0<dim*dim; g0++)
		{
			delete[] chi0s[g0];
		}
		delete[] chi0s;
		
		for(int g0=0;g0<dim*dim; g0++)
		{
			delete[] chi0[g0];
		}
		delete[] chi0;
		
		for(int g0=0;g0<dim; g0++)
		{
			delete[] chi0_gg[g0];			 			 
		}
		delete[] chi0_gg;
		
		for(int i=0;i<nomega; i++)
		{
			delete[] T[i];
		}
		delete[] T;		
		
		for(int g0=0;g0<dim; g0++)
		{
			delete[] rpa[g0];
		}
		delete[] rpa;
		
		if(system == "surface")
		{
			for(int g0=0; g0<dim_para; g0++)
			{
				delete[] chi0_para[g0];
			}
			delete[] chi0_para;
			
			for(int g0=0; g0<dim_para; g0++)
			{
				delete[] kernel[g0];
			}
			delete[] kernel;
			
			for(int g0=0; g0<dim_para; g0++)
			{
				delete[] rpa1[g0];
			}
			delete[] rpa1;  
			
			for(int g0=0; g0<dim_para; g0++)
			{
				delete[] chi_para[g0];
			}
			delete[] chi_para;
			
			delete[] flag1;
			
			for(int g0=0; g0<dim; g0++)
			{
				delete[] para_g[g0];
			}
			delete[] para_g;		
 		}
		else
		{
			for(int g0=0;g0<dim; g0++)
			{
				delete[] kernel[g0];
			}	
			delete[] kernel;
			
			for(int g0=0; g0<dim_para; g0++)
			{
				delete[] chi[g0];
			}
			delete[] chi;
		}
                				
		delete[] num_G_core;
		delete[] num_G_dis;
		delete[] G_r_core;
		delete[] num_Gvector_core;
		delete[] num_Gvector_dis;
		delete[] G_r;
		delete[] Gvec_core;
		delete[] Gvec;
		delete[] all_gcar;
		delete[] flag;
		
		cout <<"delete OK"<<endl;		
	}
	
	return;
}

//----------------------------------------------------------
//    parallel in G     collection and sorting
//----------------------------------------------------------
void Chi0_hilbert::Parallel_G()
{
	TITLE("Chi0_hilbert","Parallel_G");
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
		Gvec_core[3*g0] = pw.get_G_cartesian_projection(g0 , 0);
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

void Chi0_hilbert::Cal_Psi(int iq, complex<double> **psi_r)
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


void Chi0_hilbert::Cal_Psi_down(int iq, complex<double> **psi_r)
{
	double phase_x, phase_xy, phase_xyz;
	complex<double> exp_tmp;
	for(int ib = 0; ib < NBANDS; ib++)
	{
		ZEROS( UFFT.porter, (pw.nrxx) );
		for(int ig = wf.npwx; ig < wf.npwx + kv.ngk[iq] ; ig++)
		{
			UFFT.porter[ pw.ig2fftw[wf.igk(iq,ig - wf.npwx)] ] = wf.evc[iq](ib,ig);
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

//---------------------------------------------------------------------------
//  calculation in lcao(about to do)
//---------------------------------------------------------------------------

void Chi0_hilbert::Cal_lcao_psi()
{
}

//---------------------------------------------------------------------------
//  calculate the b matrix
//---------------------------------------------------------------------------
void Chi0_hilbert::Cal_b_lcao(int iq, int ik, int iqk)
{
	TITLE("Chi0_hilbert","Cal_b_lcao");
	double arg;
	Vector3<double> qk;
	qk = kv.kvec_c[iq] + kv.kvec_c[ik];
	//cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<endl;
	
	for(int ig=0; ig<dim; ig++)
	{
		for(int ib1=0; ib1<oband; ib1++)
		{
			ZEROS( b[ig][ib1], NBANDS);
		}
	}
	
	vector<vector<vector<complex<double>>>> phase(NLOCAL,
		vector<vector<complex<double>>>(NLOCAL,
			vector<complex<double>>(NR)));	// Peize Lin change ptr to vector at 2020.01.31
	
	int MR = 0;
	for(int iw1=0; iw1<NLOCAL; iw1++)
	{
		for(int iw2=0; iw2<NLOCAL; iw2++)
		{
			for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
			{
				arg = qk * Rcar[iw1][iw2][ir] * TWO_PI;
				phase[iw1][iw2][ir] = complex<double>( cos(arg),  sin(arg) );
				MR++;
				//cout << "phase["<<iw1<<"]["<<iw2<<"]["<<ir<<"] = "<<phase[iw1][iw2][ir]<<endl;
			}
		}
	}
	
	//------------------------------------------------
	// formula
	//------------------------------------------------
	/*for(int ig=0; ig<dim; ig++)
	{
		for(int ib1=0; ib1< oband; ib1++)
		{
			for(int ib2=0; ib2<NBANDS; ib2++) 
			{
				//------------------------------
				// calculate the matrix
				//------------------------------
				for(int iw1=0; iw1<NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<NLOCAL; iw2++)
					{
						for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
						{
							b[ig][ib1][ib2] += conj(LOWF.WFC_K[ik][ib1][iw1]) * LOWF.WFC_K[iqk][ib2][iw2] * overlap[iw1][iw2][ig][ir] * phase[iw1][iw2][ir];
						}
					}
				}				 					
			}
		}
	}*/

	vector<vector<complex<double>>> left(oband,
		vector<complex<double>>(MR));				// Peize Lin change ptr to vector at 2020.01.31
	vector<vector<complex<double>>> right(MR,
		vector<complex<double>>(NBANDS));			// Peize Lin change ptr to vector at 2020.01.31
	vector<vector<complex<double>>> b_single(oband,
		vector<complex<double>>(NBANDS)); 			// Peize Lin change ptr to vector at 2020.01.31
	for(int ig=0; ig<dim; ig++)
	{
		for(int ib1=0; ib1<oband; ib1++)
		{
			int r=0;
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
					{
						left[ib1][r] = conj(LOWF.WFC_K[ik][ib1][iw1]) * overlap[iw1][iw2][ig][ir] * phase[iw1][iw2][ir]; 
						r++;
					}
				}
			}
		}
		
		for(int ib2=0; ib2<NBANDS; ib2++)
		{
			int r=0;
			for(int iw1=0; iw1<NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<NLOCAL; iw2++)
				{
					for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
					{
						right[r][ib2] = LOWF.WFC_K[iqk][ib2][iw2];
						r++;
					}
				}
			}
		}
		
		complex<double> alpha=1;
		complex<double> beta=0;
		char transa = 'N';
		char transb = 'N';
		int M = NBANDS;
		int N = oband;
		int K = MR;
		int lda =NBANDS;
		int ldb = MR;
		int ldc =NBANDS;
		zgemm_(&transa, &transb, &M, &N, &K, &alpha, VECTOR_TO_PTR(right[0]), &lda, 
			VECTOR_TO_PTR(left[0]), &ldb, &beta,  VECTOR_TO_PTR(b_single[0]), &ldc);
		
		for(int ib1=0; ib1<oband; ib1++)
		{
			for(int ib2=0; ib2<NBANDS; ib2++)
			{
				b[ig][ib1][ib2] = b_single[ib1][ib2];
			}
		}
	}
	
	return;
}

void Chi0_hilbert::Cal_b(int iq, int ik, int iqk, int ispin)
{
	TITLE("Chi0_hilbert","Cal_b");
	Vector3<double> qk;
	qk = kv.kvec_c[iq] + kv.kvec_c[ik];
	//cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<endl;
	double phase_x, phase_xy, phase_xyz;
	Vector3<double> q = kv.kvec_d[iq];
	complex<double> exp_tmp;
	
	if(ispin == 0)
	{
		Cal_Psi(ik, psi_r1);
		Cal_Psi(iqk, psi_r2);
	}
	else if(ispin == 1)
	{
		Cal_Psi_down(ik, psi_r1);
		Cal_Psi_down(iqk, psi_r2);		
	}

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
				b[g0][ib1][ib2] = b_order[g0];
			}

		}
	}
	
	//---------------------------------------------------------------------------
	//  the gamma point can be calculated in epsilon.cpp not here!!!!!!!!!!!!!
	//---------------------------------------------------------------------------
	//if(iq == 0)
	//{
	//	for(int ib1=0; ib1<oband; ib1++)
	//	{
	//		for(int ib2=0; ib2<NBANDS; ib2++)
	//		{
	//			b[i][ib1][ib2] = complex<double>(0.0,0.0);
	//		}
	//	}
		
	//	for(int ib1=0; ib1<oband; ib1++)
	//	{
	//		for(int ib2=0; ib2<NBANDS; ib2++)
	//		{
	//			for(int ig =0; ig< kv.ngk[ik]; ig++)
	//			{
	//				b[ig][ib1][ib2] += conj(wf.evc[ik](ib1,ig)) * sqrt(((kv.kvec_c[ik]+pw.gcar[ig])*(-TWO_PI/ucell.lat0)).norm2()) * wf.evc[ik](ib2,ig);
	//			}
	//			b[ig][ib1][ib2] /=(wf.ekb[ik][ib1] - wf.ekb[ik][ib2]);
	//		}
	//	}
	//}
	
	return;
}

//-----------------------------------------------------------------
// the basic method to calculate the b matrix(just for test  1core)
//-----------------------------------------------------------------
 
/*void Chi0_hilbert::Cal_b( int iq, int ik, int iqk)
{

        Vector3<double> qk;
        qk = kv.kvec_c[iq] + kv.kvec_c[ik];
        //cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<endl;

        double phase_x, phase_xy, phase_xyz;
        Vector3<double> q = kv.kvec_d[iq];
        cout <<"q_in = "<<q.x<<" "<<q.y<<" "<<q.z<<endl;
        complex<double> exp_tmp;
        Cal_Psi1(ik);
        Cal_Psi2(iqk);

        for(int ib1=0;ib1<oband; ib1++)
            for(int ib2=0;ib2<NBANDS; ib2++)
                for(int ig=0;ig<dim;ig++)
                {
                    b[ig][ib1][ib2].real() = 0.0;
                    b[ig][ib1][ib2].imag() = 0.0; 
                }

        //for(int ib1=0; ib1<oband; ib1++)
        //{
        //        for(int ib2=0; ib2<NBANDS; ib2++)
        //        {

        //            for(int ig=0; ig<dim; ig++)
        //            {
        //                int ir=0;
        //                for(int ix=0; ix<pw.ncx; ix++)
        //                {
        //                        phase_x = (q.x + pw.gdirect[ig].x)*ix/pw.ncx;
        //                        for(int iy=0; iy<pw.ncy; iy++)
        //                        {
        //                                phase_xy = phase_x + (q.y + pw.gdirect[ig].y)*iy/pw.ncy;
        //                                for(int iz=pw.nczp_start; iz<pw.nczp_start+pw.nczp; ++iz)
        //                                {
        //                                        phase_xyz = (phase_xy + (q.z + pw.gdirect[ig].z)*iz/pw.ncz) *TWO_PI;
        //                                        exp_tmp.real() = cos(-phase_xyz);
        //                                        exp_tmp.imag() = sin(-phase_xyz);
        //                                        b[ig][ib1][ib2] = b[ig][ib1][ib2] +  conj(psi_r1[ib1][ir]) * psi_r2[ib2][ir] *exp_tmp;
        //                                        ++ir;
        //                                }
        //                        }
        //                }

        //           }
        //        }
        //}


        for(int ig=0;ig<dim;ig++)
        {
            int ir=0;
            for(int ix=0;ix<pw.ncx;ix++)
            {
                phase_x = (q.x + pw.gdirect[ig].x)*ix/pw.ncx;
                for(int iy=0; iy<pw.ncy; iy++)
                {
                    phase_xy = phase_x + (q.y + pw.gdirect[ig].y)*iy/pw.ncy;
                    for(int iz=pw.nczp_start; iz<pw.nczp_start+pw.nczp; ++iz)
                    {
                        phase_xyz = (phase_xy + (q.z + pw.gdirect[ig].z)*iz/pw.ncz) *TWO_PI;
                        exp_tmp.real() = cos(-phase_xyz);
                        exp_tmp.imag() = sin(-phase_xyz);
                        for(int ib1=0;ib1<oband;ib1++)
                        {
                            for(int ib2=0;ib2<NBANDS;ib2++)
                            {
                                b[ig][ib1][ib2] = b[ig][ib1][ib2] +  conj(psi_r1[ib1][ir]) * psi_r2[ib2][ir] *exp_tmp;
                            }
                        }
                        ++ir;
                    }
                }
            }
        }


        for(int ib1=0;ib1<oband; ib1++)
            for(int ib2=0;ib2<NBANDS; ib2++)
                for(int ig=0;ig<dim;ig++)
                {
                    b[ig][ib1][ib2] = b[ig][ib1][ib2]/pw.ncxyz; 
                }

        return;

}*/

// calculate the chi0s matrix
void Chi0_hilbert::Cal_Chi0s(int iq)
{
	TITLE("Chi0_hilbert","Cal_Chi0s");
	double delta_e,e1,e2;
	complex<double> weight1, weight2;
	for(int g0=0; g0<dim; g0++)
	{
		for(int g1=0; g1<dim; g1++)
		{
			ZEROS( chi0s[g1+g0*dim], nomega);
		}
	}
	
	//---------------------------------------------------------------------------------
	// some test(use our program's wavefunction and other program's occupation to check)
	//---------------------------------------------------------------------------------
	/*double Q[kv.nks][NBANDS];
	ifstream ifs("band.dat");
	if(!ifs)
	{
		WARNING_QUIT("chi0_hilbert","Can't find band.dat");
	}
	for(int i=0;i<kv.nks;i++)
	{
		for(int j=0;j<NBANDS;j++)
		{
			ifs >> Q[i][j];
		}
	}*/
	
	//---------------------------------------------------------------------------------
	// no spin
	//---------------------------------------------------------------------------------
	if(NSPIN == 1) 
	{
		for(int ik=0; ik<kv.nks; ik++)
		{
			int iqk = Cal_iq(ik, iq, kv.nmp[0], kv.nmp[1], kv.nmp[2]);
			if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
			{
				Cal_b(iq, ik, iqk, 0);
			}
			else
			{
				Cal_b_lcao(iq, ik, iqk);
			}

			for(int ib1=0; ib1<oband; ib1++)
			{
				for(int ib2=0; ib2<NBANDS; ib2++)
				{
					delta_e = wf.ekb[iqk][ib2] - wf.ekb[ik][ib1];
					//delta_e = Q[ik][ib2] - Q[ik][ib1];
					if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						weight2 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/ucell.omega, 0.0);
						for(int g0=0; g0<dim; g0++)
						{
							for(int g1=0; g1<dim; g1++)
							{
								chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
								chi0s[g1+g0*dim][n+1] += weight2 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
							}
						}
					}
					else if((delta_e > ((nomega-1) * domega) || delta_e == ((nomega-1) * domega)) && delta_e < (nomega * domega))
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						for(int g0=0; g0<dim; g0++)
						{
							for(int g1=0; g1<dim; g1++)
							{
								chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
							}
						}
					}
				}
			}
		}
	}
	//---------------------------------------------------------------------------------
	// spin = 2 ( fermi_level != 0 haven't been tested)
	//---------------------------------------------------------------------------------
	else if(NSPIN == 2)
	{
		// spin up
		for(int ik=0; ik<kv.nks/2; ik++)
		{
			int iqk = Cal_iq(ik, iq, kv.nmp[0], kv.nmp[1], kv.nmp[2]);
			if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
			{
				Cal_b(iq, ik, iqk, 0);
			}
			else
			{
				Cal_b_lcao(iq, ik, iqk);
			}
			
			for(int ib1=0; ib1<oband; ib1++)
			{
				for(int ib2=0; ib2<NBANDS; ib2++)
				{
					delta_e = wf.ekb[iqk][ib2] - wf.ekb[ik][ib1];
					if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						//weight1 = complex<double>( (wf.wg(ik,ib1) - wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/ucell.omega, 0.0);                        
						//weight2 = complex<double>( (wf.wg(ik,ib1) - wf.wg(iqk,ib2)) * (delta_e - e1)/domega/ucell.omega, 0.0);
						weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						weight2 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/ucell.omega, 0.0);
						for(int g0=0; g0<dim; g0++)
						{
							for(int g1=0; g1<dim; g1++)
							{
								chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
								chi0s[g1+g0*dim][n+1] += weight2 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
							}
						}
					}
					else if((delta_e > ((nomega-1) * domega) || delta_e == ((nomega-1) * domega)) && delta_e < (nomega * domega))
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						//weight1 = complex<double>( (wf.wg(ik,ib1) - wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						for(int g0=0; g0<dim; g0++)
						{
							for(int g1=0; g1<dim; g1++)
							{
								chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
							}
						}
					}
				}
			}			
		}
		// spin down
		for(int ik=kv.nks/2 + 1; ik<kv.nks; ik++)
		{
			int iqk = Cal_iq(ik-kv.nks/2, iq, kv.nmp[0], kv.nmp[1], kv.nmp[2]) + kv.nks/2;
			if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
			{
				Cal_b(iq+kv.nks/2, ik, iqk, 0);
			}
			else
			{
				Cal_b_lcao(iq+kv.nks/2, ik, iqk);
			}
			for(int ib1=0; ib1<oband; ib1++)
			{
				for(int ib2=0; ib2<NBANDS; ib2++)
				{
					delta_e = wf.ekb[iqk][ib2] - wf.ekb[ik][ib1];
					if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						//weight1 = complex<double>( (wf.wg(ik,ib1) - wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						//weight2 = complex<double>( (wf.wg(ik,ib1) - wf.wg(iqk,ib2)) * (delta_e - e1)/domega/ucell.omega, 0.0);  
						weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						weight2 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/ucell.omega, 0.0);						
						for(int g0=0; g0<dim; g0++)
						{
							for(int g1=0; g1<dim; g1++)
							{
								chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
								chi0s[g1+g0*dim][n+1] += weight2 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
							}
						}
					}
					else if((delta_e > ((nomega-1) * domega) || delta_e == ((nomega-1) * domega)) && delta_e < (nomega * domega))
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						//weight1 = complex<double>( (wf.wg(ik,ib1) - wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
						for(int g0=0; g0<dim; g0++)
						{
							for(int g1=0; g1<dim; g1++)
							{
								chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
							}
						}
					}
				}
			}			
		}
	}
	else if(NSPIN == 4)
	{
		cout<<"NSPIN = "<<NSPIN<<endl;
		for(int ik=0; ik<kv.nks; ik++)
		{
			int iqk = Cal_iq(ik, iq, kv.nmp[0], kv.nmp[1], kv.nmp[2]);
			for(int ispin =0; ispin<2; ispin++)
			{
				if(BASIS_TYPE == "pw" || BASIS_TYPE == "lcao_in_pw")
				{
					Cal_b(iq, ik, iqk, ispin);
				}
				else
				{
					Cal_b_lcao(iq, ik, iqk);
				}
				cout<<"ik = "<<ik<<" ispin = "<<ispin<<endl;
				for(int g=0; g<dim; g++)
				{
					cout<<"b["<<g<<"][2][3]"<<" "<<b[g][2][3]<<endl;
				}	

				for(int ib1=0; ib1<oband; ib1++)	
				{
					for(int ib2=0; ib2<NBANDS; ib2++)
					{
						delta_e = wf.ekb[iqk][ib2] - wf.ekb[ik][ib1];
						//delta_e = Q[ik][ib2] - Q[ik][ib1];
						if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
						{
							int n = int(delta_e/domega);
							e1 = double(n) * domega;
							e2 = double(n+1) * domega;
							weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
							weight2 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/ucell.omega, 0.0);
							for(int g0=0; g0<dim; g0++)
							{
								for(int g1=0; g1<dim; g1++)
								{
									chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
									chi0s[g1+g0*dim][n+1] += weight2 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
								}
							}
						}
						else if((delta_e > ((nomega-1) * domega) || delta_e == ((nomega-1) * domega)) && delta_e < (nomega * domega))
						{
							int n = int(delta_e/domega);
							e1 = double(n) * domega;
							e2 = double(n+1) * domega;
							weight1 = complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/ucell.omega, 0.0);
							for(int g0=0; g0<dim; g0++)
							{
								for(int g1=0; g1<dim; g1++)
								{
									chi0s[g1+g0*dim][n] += weight1 * b[g0][ib1][ib2] * conj(b[g1][ib1][ib2]);
								}
							}
						}
					}
				}
			}
		}
	}
	
	return;
}

// calculte the transportation matrix
void Chi0_hilbert::Cal_T()
{
	TITLE("Chi0_hilbert","Cal_T");
	complex<double> M1, M2;
	for(int n1=0; n1<nomega; n1++)
	{
		for(int n=0; n<nomega; n++)
		{
			double n_e = double(n) * domega;
			double n1_e = double(n1) * domega;
			M1 = complex<double>( (n_e - n1_e), eta);
			M2 = complex<double>( (n_e + n1_e), eta);
			T[n1][n] = 1.0/M1 - 1.0/M2;
		}
	}
	
	return;
}

// calculate the chi0
void Chi0_hilbert::Cal_Chi0()
{
	TITLE("Chi0_hilbert","Cal_Chi0");
	vector<vector<complex<double>>> chi0s2(dim*dim,
		vector<complex<double>>(nomega));			// Peize Lin change ptr to vector at 2020.01.31
	vector<vector<complex<double>>> T2(nomega,
		vector<complex<double>>(nomega));			// Peize Lin change ptr to vector at 2020.01.31
	vector<vector<complex<double>>> chi02(dim*dim,
		vector<complex<double>>(nomega));			// Peize Lin change ptr to vector at 2020.01.31
	for(int i=0;i<dim*dim;i++)
	{
		for(int j=0;j<nomega;j++)
		{
			chi0s2[i][j] = chi0s[i][j];
		}
	}
	
	for(int i=0;i<nomega;i++)
	{
		for(int j=0;j<nomega;j++)
		{
			T2[i][j] = T[i][j];
		}
	}
	
	//---------------------------------------------------------------------------------
	// use lapack to reduce the time
	//---------------------------------------------------------------------------------
	complex<double> alpha=1;
	complex<double> beta=0;
	char transa = 'N';
	char transb = 'N';
	int M = nomega;
	int N = dim*dim;
	int K = nomega;
	int lda =nomega;
	int ldb =nomega;
	int ldc =nomega;
	
	zgemm_(&transa, &transb, &M, &N, &K, &alpha, VECTOR_TO_PTR(T2[0]), &lda, VECTOR_TO_PTR(chi0s2[0]), &ldb, &beta, VECTOR_TO_PTR(chi02[0]), &ldc);
	
	for(int i=0;i<dim*dim;i++)
	{
		for(int j=0;j<nomega;j++)
		{
			chi0[i][j] = chi02[i][j];
		}
	}
	
	/*for(int i=0;i<dim*dim;i++)
	{
		ZEROS(chi0[i], nomega);
		for(int j=0;j<nomega;j++)
		{
			for(int k=0;k<nomega;k++)
			{
				chi0[i][j] = chi0[i][j] + chi0s[i][k] * f(k,j);
			}
		}
	}*/
	
	return;
}

complex<double> Chi0_hilbert::f(int k, int j)
{
	complex<double> m,n,l;
	m = complex<double>( (j - k) * domega, eta);
	n = complex<double>( (j + k) * domega, eta);
	l = 1.0/m - 1.0/n;
	return l; 	            
}


// calculate the chi
void Chi0_hilbert::Cal_Chi(int iq)
{
	TITLE("Chi0_hilbert","Cal_Chi");
	double q = sqrt(((kv.kvec_c[iq])*(TWO_PI/ucell.lat0)).norm2());
	vector<complex<double>> gather_chi(nomega);			// Peize Lin change ptr to vector at 2020.01.31
	for(int i=0; i<nomega; i++)
	{
		for(int g0=0; g0<dim; g0++)
		{
			for(int g1=0; g1<dim; g1++)
			{
				chi0_gg[g0][g1] = chi0[g1+g0*dim][i];
			}
		}
		
		Cal_Rpa(iq);
		//cout <<"rpa ok"<<endl;
		
		for(int g0=0; g0<dim; g0++)
		{
			ZEROS( chi[g0], dim);
		}
		
		//---------------------------------------------------------------------------------
		//  an approximate method:  Now we calculate epsilon^-1 = 1 + v * chi = 
		//  1+ v*[ chi0 * (1 - v * chi0)^-1] ~= (1-v*chi0)^-1 
		//  I tested it, and the two method gives the similar results.
		//---------------------------------------------------------------------------------
		/*for(int g0 =0; g0<dim; g0++)
		{
			for(int g1 =0; g1<dim; g1++)
			{
				chi[g0][g1] = complex<double>( rpa[g0][g1].real(), rpa[g0][g1].imag());
			}
		}
		cout <<"iq = "<<iq<< " epsilon["<<i<<"] = "<<chi[0][0].real()<<"  "<<chi[0][0].imag() <<endl;*/
		
		CMatrixMul(dim, dim, dim, chi0_gg, rpa, chi);
		//cout <<"chi ok"<<endl;
		
		gather_chi[i] = chi[0][0];
		
		cout <<"iq = "<<iq<< " chi["<<i<<"] = "<<chi[0][0].real()<<"  "<<chi[0][0].imag() <<endl;
		//cout << "coulomb_cutoff = " << coulomb_cutoff << endl;
		//---------------------------------------------------------------------
		// -Im{epsilon^-1} = -4pi/|q|^2 * Im{chi}
		//---------------------------------------------------------------------
		if(!coulomb_cutoff)
		{
			cout <<"iq = "<<iq<<" epsilon["<<i<<"] = "<<8*PI/q/q*chi[0][0].real()<<"  "<<8*PI/q/q*chi[0][0].imag()<<endl;
		}
		else
		{
			//---------------------------------------------------------------------
			// -Im{epsilon^-1} = -4pi/|q|^2 * (1 - exp^{-q*lz/2}) * Im{chi}
			//---------------------------------------------------------------------
			double factor;
			factor = 1.0 - exp(-q * (ucell.latvec.e33 * ucell.lat0)/2.0);
			cout <<"iq = "<<iq<<" epsilon["<<i<<"] = "<<8*PI/q/q * factor * chi[0][0].real()<<"  "<<8*PI/q/q * factor * chi[0][0].imag()<<endl; 
		}
		
		//--------------------------------------------------------------------
		//  Re{epsilon}  NJP(2014)
		//--------------------------------------------------------------------
		//complex<double> M;
		//M = 1.0/(8*PI/q/q*chi[0][0]);
		//cout <<"iq = "<<iq<<" epsMreal["<<i<<"] = "<<M.real()<<endl;
		//cout <<"iq = "<<iq<<" epsMimag["<<i<<"] = "<<M.imag()<<endl; 
	}
	
	if(out_chi)
	{
		stringstream sschi;
		sschi << global_out_dir << "chi_"<<iq<<".dat";
		ofstream ofschi(sschi.str().c_str());
		for(int i=0;i<nomega;i++)
		{
			ofschi << i * domega <<"   "<<gather_chi[i].real()<<"   "<<gather_chi[i].imag() <<endl;
		}
		ofschi.close();
	}
	
	if(out_epsilon)
	{
		stringstream sseps;
		sseps << global_out_dir << "Imeps^-1_"<<iq<<".dat";
		ofstream ofseps(sseps.str().c_str());
		ofseps<<"Energy(Ry)"<<"   "<<"-Im{epsilon^-1}"<<endl;
		if(!coulomb_cutoff)
		{
			for(int i=0;i<nomega;i++)
			{
				ofseps << i * domega <<"   "<<-8*PI/q/q*gather_chi[i].imag()<<endl;
			}
		}
		else
		{
			double factor;
			factor = 1.0 - exp(-q * (ucell.latvec.e33 * ucell.lat0)/2.0);
			for(int i=0;i<nomega;i++)
			{
				ofseps << i * domega <<"   "<<-8*PI/q/q*factor*gather_chi[i].imag()<<endl;
			}			
		}
		ofseps.close();
	}
	
	return;
}

void Chi0_hilbert::Cal_Chi_surface(int iq)
{
	TITLE("Chi0_hilbert","Cal_Chi_surface");
	complex<double> g;
	for(int i=0; i<nomega; i++)
	{
		for(int g0=0; g0<dim; g0++)
		{
			for(int g1=0; g1<dim; g1++)
			{
				chi0_gg[g0][g1] = chi0[g1+g0*dim][i];
			}
		}
		
		chi0_para_g();
		cout <<"iq = "<<iq<< " chi0_para["<<i<<"] = "<<chi0_para[0][0].real()<<"  "<<chi0_para[0][0].imag() <<endl;
		Cal_Rpa1();
		cout <<"iq = "<<iq<< " epsilon["<<i<<"] = "<<rpa1[0][0].real()<<"  "<<rpa1[0][0].imag() <<endl;
		
		//g = Cal_g(iq);
		//cout<<"iq = "<<iq<<" g["<<i<<"] = "<<g.real()<<"  "<<g.imag()<<endl;
	}
	
	return;
}

// calculate kernel 
void Chi0_hilbert:: Cal_kernel(int iq)
{
	if(kernel_type == "rpa")
	{
		//-------------------------------------------------------------
		//   Coulomb   4PI/|q+G|^2
		//-------------------------------------------------------------
		if(!coulomb_cutoff)
		{
			for(int g0=0; g0<dim; g0++)
				for(int g1=0; g1<dim; g1++)
				{
					if(g0 == g1)
					{
						kernel[g0][g1] = 8.0 * PI/qg2(iq,g0);
					}
					else
					{
						kernel[g0][g1] = complex<double>(0.0, 0.0);
					}
				}
		}
		//-------------------------------------------------------------
		//  Coulomb cutoff: use in surface like Graphene
		//-------------------------------------------------------------
		else
		{
			double sign;
			for(int g0=0; g0<dim; g0++)
			{
				int nz = int(fabs(all_gcar[g0].z * ucell.latvec.e33) + 0.5 );
				sign = (nz%2==0)?1.0:-1.0;
				for(int g1=0; g1<dim; g1++)
				{
					if(g0 == g1)
					{
						kernel[g0][g1] = 8.0 * PI/qg2(iq,g0) * (1 - sign * exp(-qG(iq,g0)*(ucell.latvec.e33 * ucell.lat0)/2.0));
					}
					else
					{
						kernel[g0][g1] = complex<double>(0.0, 0.0);
					}
				}
			}			
		}
	}
	else if(kernel_type == "tdlda")
	{
		
	}
}
// calculate (1 -vx) matrix
void Chi0_hilbert::Cal_Rpa(int iq)
{
	//double rc3, rc, factor, qg;
	//-------------------------------------------------------------
	//  Coulomb   4PI/|q+G|^2
	//-------------------------------------------------------------
	/*if(!coulomb_cutoff)
	{
		for(int g0=0; g0<dim; g0++)
		{
			for(int g1=0; g1<dim; g1++)
			{
				if (g0 != g1)
				{
					rpa[g0][g1] = -8.0 * PI/qg2(iq,g0) * chi0_gg[g0][g1];
				}
				else
				{
					rpa[g0][g1] = 1.0 - 8.0 * PI/qg2(iq,g0) * chi0_gg[g0][g1];
				}
			}
		}
	}*/
	//-------------------------------------------------------------
	//  Coulomb cutoff: use in surface like Graphene
	//-------------------------------------------------------------
	/*else
	{
		double sign;
		for(int g0=0; g0<dim; g0++)
		{
			int nz = int(fabs(all_gcar[g0].z * ucell.latvec.e33) + 0.5 );
			sign = (nz%2==0)?1.0:-1.0;
			//cout <<"g0= "<<g0<<" nz= "<<nz<<"  sign= "<<sign<<endl;
			for(int g1=0; g1<dim; g1++)
			{
				if(g0 != g1)
				{
					rpa[g0][g1] = -8.0 * PI/qg2(iq,g0) * (1 - sign * exp(-qG(iq,g0)*(ucell.latvec.e33 * ucell.lat0)/2.0)) * chi0_gg[g0][g1];
				}
				else
				{
					rpa[g0][g1] = 1.0 - 8.0 * PI/qg2(iq,g0) * (1 - sign * exp(-qG(iq,g0)*(ucell.latvec.e33 * ucell.lat0)/2.0)) * chi0_gg[g0][g1];
				}
			}
			
		}
	}*/

	CMatrixMul(dim, dim, dim, kernel, chi0_gg, rpa);
	for(int g0=0; g0<dim; g0++)
		for(int g1=0; g1<dim; g1++)
		{
			if(g0 == g1)
			{
				rpa[g0][g1] = 1.0 - rpa[g0][g1];
			}
			else
			{
				rpa[g0][g1] = -rpa[g0][g1];
			}
		}
	
	//-------------------------------------------------------------
	//  a method used in PRB 90,161410
	//  they calculate epsilon not epsilon^-1 to get the peak 
	//  here we use lapack to get the eigenvalues
	//-------------------------------------------------------------
	/*complex<double> a[dim*dim];
	for(int g0=0;g0<dim;g0++)
	{
		for(int g1=0;g1<dim;g1++)
		{
			for(int g1=0;g1<dim;g1++)
		}
	}
	
	for(int g0=0;g0<dim*dim;g0++)
	{
		cout << a[g0]<<"  ";
	}
	cout <<endl;

	char jobvl = 'V';
	char jobvr = 'V';
	int N = dim;
	int lda = dim;
	complex<double>* w = (complex<double>*)malloc( sizeof(complex<double>) * N);
	int ldvr = dim;
	complex<double>* vr = (complex<double>*)malloc( sizeof(complex<double>) * N * ldvr);
	int ldvl = dim;
	complex<double>* vl = (complex<double>*)malloc( sizeof(complex<double>) * N * ldvl);
	int lwork = N * 4;
	double *rwork = (double*)malloc( sizeof(double) * 2 * N);
	complex<double> *work = (complex<double>*)malloc( sizeof(complex<double>) * lwork);
	int info;

	zgeev_(&jobvl, &jobvr, &N, a, &lda, w, vl , &ldvl , vr, &ldvr, work, &lwork, rwork, &info);

	cout <<"eigen  ";
	for(int g0=0;g0<dim;g0++)
	{
		cout<<w[g0].real()<<"   "<<w[g0].imag()<<"   ";
	}
	cout<<endl;*/
	
	//-------------------------------------------------------------
	// inverse the (I -vx)
	//-------------------------------------------------------------
	int l = Cinv(dim, rpa);
	if (l == 0)
	{
		WARNING_QUIT("chi0_hilbert","(I-v*chi0) is a singular matrix !!");
	}
	
	//-------------------------------------------------------------
	//  inverse the matrix with lapack, but it can only deal with 
	//  the matrix which dimension is under 4, some bugs exist!!!!
	//-------------------------------------------------------------
	//C_inverse(dim, rpa);
	
	return;
}

double Chi0_hilbert::qg2( int iq, int g0)
{
	double qg2;
	qg2 = ((kv.kvec_c[iq]+all_gcar[g0])*(TWO_PI/ucell.lat0)).norm2();
	
	return qg2;
}

//-------------------------------------------------------------
// inverse a complex matrix by Gauss elimination method
//-------------------------------------------------------------
int Chi0_hilbert::Cinv(int n, complex<double>** a)
{
	int i,j,k;	//l,u,v,w;
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
		
		//l=k*n+k;
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

//-------------------------------------------------------------
// Matrix multiplication 
//-------------------------------------------------------------
void Chi0_hilbert::CMatrixMul(int m, int n, int l, complex<double>** A, complex<double>** B, complex<double>** C)
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

//-------------------------------------------------------------
//  get the iqk(= iq + ik) in uniform kpoints
//-------------------------------------------------------------

int Chi0_hilbert::Cal_iq(int ik, int iq, int a, int b, int c)
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

//-------------------------------------------------------------
//  calculate g            exist bugs!!!!!!!!!!!
//-------------------------------------------------------------
int Chi0_hilbert::parallel_g()
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

void Chi0_hilbert::chi0_para_g()
{
	for(int g0=0; g0<dim_para; g0++)
	{
		ZEROS( chi0_para[g0], dim_para);
	}
	
	for(int g0=0; g0<dim; g0++)
	{
		for(int g1=0; g1<dim; g1++)
		{
			chi0_para[flag1[g0]][flag1[g1]] += chi0_gg[g0][g1];
		}
	}
	
	return;
}

void Chi0_hilbert::Cal_kernel_2D(int iq)
{
	if(kernel_type == "rpa")
	{
		for(int g0=0; g0<dim_para; g0++)
		{
			ZEROS( kernel[g0], dim_para);
		}
		
		for(int g0=0;g0<dim_para;g0++)
			for(int g1=0;g1<dim_para;g1++)
			{
				if(g0 == g1)
				{
					kernel[g0][g1] = 4* PI/qg(iq,g0) * (ucell.latvec.e33 * ucell.lat0);
				}
				else
				{
					kernel[g0][g1] = complex<double>(0.0,0.0);
				}
			}
	}
	else if(kernel_type == "tdlda")
	{
		
	}
	
	return;
}

double Chi0_hilbert::qG(int iq, int g0)
{
	double qG;
	qG = sqrt((kv.kvec_c[iq].x + all_gcar[g0].x) * (kv.kvec_c[iq].x + all_gcar[g0].x) + (kv.kvec_c[iq].y + all_gcar[g0].y) * (kv.kvec_c[iq].y + all_gcar[g0].y)) * (TWO_PI/ucell.lat0);
     	
	return qG;
}

double Chi0_hilbert::qg(int iq, int g0)
{
	double qg;
	qg = sqrt((kv.kvec_c[iq].x + para_g[g0][0]) * (kv.kvec_c[iq].x + para_g[g0][0]) + (kv.kvec_c[iq].y + para_g[g0][1]) * (kv.kvec_c[iq].y + para_g[g0][1])) * (TWO_PI/ucell.lat0);
	
	return qg;
}

void Chi0_hilbert::Cal_Rpa1()
{
	CMatrixMul(dim_para, dim_para, dim_para, kernel, chi0_para, rpa1); 
	for(int i=0;i<dim_para;i++)
	{
		for(int j=0;j<dim_para; j++)
		{
			if(i == j)
			{
				rpa1[i][j] = 1.0 - rpa1[i][j];
			}
			else
			{
				rpa1[i][j] = -rpa1[i][j];
			}
		}
	}
	
	int l = Cinv(dim_para, rpa1);
	if (l == 0)
	{
		WARNING_QUIT("chi0_hilbert","(I-v*chi0) is a singular matrix !!");
	}
	
	return;
}

void Chi0_hilbert::chi_para_g()
{
	CMatrixMul(dim_para, dim_para, dim_para, chi0_para, rpa1, chi_para);
	
	return;
}


complex<double> Chi0_hilbert:: Cal_g(int iq)
{	
    complex<double> g;
    double L = ucell.latvec.e33 * ucell.lat0;
    double dz = L/pw.ncz;
    double q =  sqrt(((kv.kvec_c[iq])*(TWO_PI/ucell.lat0)).norm2());
    
    g = complex<double>(0.0,0.0);

    for(int z0=0;z0<pw.ncz;z0++)
        for(int z1=0;z1<pw.ncz;z1++)
        {
            double exp_phase = exp(-q*(z0+z1)*dz);
            g += chi_para[z0][z1] * exp_phase * dz * dz;
        }
    g = -4*PI/q * g;

    return g;
}


//---------------------------------------------
//  matrix inverse( exist some bugs!)
//---------------------------------------------

/*void Chi0_hilbert:: C_inverse(int n, complex<double>** a)
{
         int info;

         //complex<double> M[n][n],N[n][n];
         complex<double> M[n*n],N[n*n];

         for(int i=0; i<n; i++)
         {
             for(int j=0; j<n; j++)
             {
                 int index = i * n + j;
                 M[index] = a[i][j];
                 //M[i][j] = a[i][j];
             }
         }
          
         cout<<"a in "<<endl;
         for(int i=0; i<n; i++)
         {
             for(int j=0; j<n; j++)
             {
                 int index = i * n + j;
                 //cout <<M[i][j]<<" ";
                 cout<<M[index]<<" ";
             }
             cout <<endl;
         }
   
         int dim = n;
         int nrhs = n;
         int lda = n;
         int ipiv = n;
         int ldb = n;
          
         //zgesv_(&dim, &nrhs, M[0], &lda, &ipiv, N[0], &ldb, &info);
         //zgesv_(&dim, &nrhs, M, &lda, &ipiv, N, &ldb, &info);
         
         zgetrf_(&dim, &dim, M, &lda, &ipiv, &info);
        
         zgetri_(&dim, M, &lda, &ipiv, N, &dim, &info);  


         cout<<"a out "<<endl;
         for(int i=0; i<n; i++)
         {
             for(int j=0; j<n; j++)
             {
                 int index = i * n + j;
                 //cout <<N[i][j]<<" ";
                 cout <<M[index]<<" ";
             }
             cout <<endl;
         }

         for(int i=0; i<n; i++)
         {
             for(int j=0; j<n; j++)
             {
                 int index = i * n + j;
                 //a[i][j] = N[i][j];
                 a[i][j] = M[index];
             }
         }
        
         //delete[] M;
         //delete[] N;
 
         return;
}*/

//----------------------------------------------------
// plot function
//----------------------------------------------------
void Chi0_hilbert::plot_chi0(int iq)
{
	stringstream ss;
	ss << global_out_dir <<"chi0"<<"_"<<iq<<".dat";
	ofstream ofs(ss.str().c_str());
	for(int i=0;i<nomega;i++)
	{
		ofs <<i*domega<<"    "<<chi0[0][i].real()<<"    "<<chi0[0][i].imag()<<endl;
	}
	ofs.close();
	
	return;
}



