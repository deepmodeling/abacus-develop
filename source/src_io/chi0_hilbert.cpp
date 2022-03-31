//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-14
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating some optical properties using the linear response
//theory . and it is a basic work and use the hilbert-transform to reduce the
//calculated quantities.
//-------------------------------------------------------------------------------
#include "../src_pw/global.h"
#include "chi0_hilbert.h"
#include "../src_pw/hamilt.h"
#include "../src_lcao/wavefunc_in_pw.h"
#include "optical.h"
#include "../src_pw/klist.h"
#include <iostream>
#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif

using namespace std;

namespace GlobalC
{
Chi0_hilbert chi0_hilbert;
}

std::complex<double> *** Chi0_hilbert::wfc_k_grid = nullptr;
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
#include "../src_pw/occupy.h"
void Chi0_hilbert::Chi()
{
	ModuleBase::TITLE("Chi0_hilbert","Chi");
	//---------------------------------------
	// check nomega
	//---------------------------------------
	double energy_scale = domega * nomega;
	std::cout <<"energy_scale = "<<energy_scale<<std::endl;
	double max_e = 0.0;
	for(int ik1=0; ik1<GlobalC::kv.nks; ik1++)
	{
		for(int ik2=0; ik2<GlobalC::kv.nks; ik2++)
		{
			if(max_e < (GlobalC::wf.ekb[ik1][GlobalV::NBANDS-1] - GlobalC::wf.ekb[ik2][0]))
			{
				max_e = GlobalC::wf.ekb[ik1][GlobalV::NBANDS-1] - GlobalC::wf.ekb[ik2][0];
			}
		}
	}
	std::cout <<"max_e = "<<max_e<<std::endl;
	if(max_e > energy_scale)
	{
		int nomega_right = int(max_e/domega);
		std::cout << "nomega must >= "<<nomega_right<<std::endl;
		ModuleBase::WARNING_QUIT("chi0_hilbert","nomega is too small!!!!!!!!!");
	}

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
	std::vector<int> qc(GlobalC::kv.nks); 			// Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<double> ql(GlobalC::kv.nks);			// Peize Lin change ptr to std::vector at 2020.01.31
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

	//----------------------------------------------------------
	//    parallel in G     collection and sorting
	//----------------------------------------------------------
	Parallel_G();

	if(system == "surface")    // calculate the g_//
	{
		dim_para = parallel_g();
		//std::cout << "dim_para = "<<dim_para<<std::endl;
	}

	//std::cout <<"GlobalC::en.ef = "<<GlobalC::en.ef<<std::endl;
	//std::cout <<"smearing_sigma = "<<INPUT.smearing_sigma<<std::endl;

	//------------------------------------------------------------------
	// change the fermi level or not
	//------------------------------------------------------------------
	double energy = GlobalC::en.ef + fermi_level;
	//std::cout <<"energy = "<<energy<<std::endl;

	cweight = new double*[GlobalC::kv.nks];
	for(int i=0;i<GlobalC::kv.nks;i++)             // new occupation
	{
		cweight[i] = new double[GlobalV::NBANDS];
	}

	for(int ik = 0;ik < GlobalC::kv.nks;ik++)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			//----------------------------------------------------------
			//  gauss smearing
			//----------------------------------------------------------
			if(INPUT.smearing_method == "gauss" || INPUT.smearing_method == "gaussian")
			{
				cweight[ik][ib] = GlobalC::kv.wk[ik] * Occupy::wgauss( (energy - GlobalC::wf.ekb[ik][ib] )/INPUT.smearing_sigma, 0);
			}
			//----------------------------------------------------------
			//  fixed smearing
			//----------------------------------------------------------
			else if(INPUT.smearing_method == "fixed")
			{
				if((GlobalC::wf.ekb[ik][ib]-energy)<-0.0001 )
				{
					cweight[ik][ib] = GlobalC::kv.wk[ik] * 1.0;
				}
				else if(fabs(GlobalC::wf.ekb[ik][ib]-energy)<0.0001 )
				{
					cweight[ik][ib] = GlobalC::kv.wk[ik] * 0.5;
				}
				else
				{
					cweight[ik][ib] = 0.0;
				}
			}
			//----------------------------------------------------------
			//   fermi-dirac smearing
			//----------------------------------------------------------
			else if(INPUT.smearing_method == "fd")
			{
				cweight[ik][ib] = GlobalC::kv.wk[ik] * Occupy::wgauss( (energy - GlobalC::wf.ekb[ik][ib] )/INPUT.smearing_sigma, -99);
			}
			else
			{
				ModuleBase::WARNING_QUIT("chi0_hilbert","Sorry! This smearing hasn't been realized here yet!!!");
			}			
		} 
	}

	//----------------------------------------------------------
	//  calculate the total electrons
	//----------------------------------------------------------
	double elec = 0.0;
	for(int ik = 0;ik < GlobalC::kv.nks;ik++)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			elec += cweight[ik][ib];
		}
	}

	std::cout <<"elec = "<<elec<<std::endl;

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
			if( cweight[ik][ib] > 0.0001)
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
		if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
		{
			for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
			{
				time_t start,end;
				start = time(NULL);
				Cal_Chi0s(iq);  // calculate chi0s
				end = time(NULL);
				std::cout<<"chi0s time:: "<<end-start<<" s"<<std::endl;

				time_t start1,end1;
				start1 = time(NULL);
				Cal_T();
				end1 = time(NULL);
				std::cout<<"T time:: "<<end1-start1<<" s"<<std::endl;

				time_t start2,end2;
				start2 = time(NULL);
				Cal_Chi0();   // calculate  chi0
				end2 = time(NULL);
				std::cout<<"chi0 time:: "<<end2-start2<<" s"<<std::endl;

				if(out_chi0)
				{
					plot_chi0(iq);
				}
				for(int i=0; i<nomega; i++)
				{
					std::cout <<"iq = "<<iq<<" chi0["<<i<<"] = "<<chi0[0][i].real()<<"  "<<chi0[0][i].imag()<<std::endl;
				}

				Cal_kernel(iq);

				time_t start3,end3;
				start3 = time(NULL);
				Cal_Chi(iq);  // calculate chi
				end3 = time(NULL);
				std::cout<<"chi time:: "<<end3-start3<<" s"<<std::endl;
			}
		}
		//-------------------------------------------------
		// lcao method
		//-------------------------------------------------
		else if(GlobalV::BASIS_TYPE == "lcao")
		{
			//-------------------------------------------------
			// first, readin files( "nearest.dat" and "overlap_q")
			//-------------------------------------------------
			int tmp; int icount = 0;
			//------------------------
			// nearest.dat
			//------------------------
			std::stringstream ss;
			ss << GlobalV::global_out_dir <<"nearest.dat";
			std::cout << ss.str().c_str() << std::endl;
			std::ifstream ifsn(ss.str().c_str());
			if(!ifsn)
			{
				ModuleBase::WARNING_QUIT("chi0_hilbert", "Can't find the nearest.dat file!");
			}
			ifsn >> NR;
			std::cout << "NR in = "<<NR << std::endl;
			
			R = new ModuleBase::Vector3<int>** [GlobalV::NLOCAL];
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				R[iw1] = new ModuleBase::Vector3<int>* [GlobalV::NLOCAL];
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					R[iw1][iw2] = new ModuleBase::Vector3<int>[NR];
				}
			}
			
			Rcar = new ModuleBase::Vector3<double>** [GlobalV::NLOCAL];
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				Rcar[iw1] = new ModuleBase::Vector3<double>* [GlobalV::NLOCAL];
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					Rcar[iw1][iw2] = new ModuleBase::Vector3<double>[NR];
				}
			}

			Rmax = new int* [GlobalV::NLOCAL];
			for(int iw=0; iw<GlobalV::NLOCAL; iw++)
			{
				Rmax[iw] = new int[GlobalV::NLOCAL];
			}

			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					//std::cout <<"iw1 = "<<iw1<<" iw2 = "<<iw2<<std::endl;
					ifsn >> tmp;
					ifsn >> tmp;
					ifsn >> Rmax[iw1][iw2]; //std::cout << Rmax[iw1][iw2] << std::endl;
					for(int i=0; i<Rmax[iw1][iw2]; i++)
					{
						ifsn >> R[iw1][iw2][i].x; //std::cout <<R[iw1][iw2][i].x;
						ifsn >> R[iw1][iw2][i].y; //std::cout <<R[iw1][iw2][i].y;
						ifsn >> R[iw1][iw2][i].z; //std::cout <<R[iw1][iw2][i].z; std::cout << std::endl;
					}
				}
			}
			ifsn.close();

			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					for(int i=0; i<Rmax[iw1][iw2]; i++)
					{
						Rcar[iw1][iw2][i].x = GlobalC::ucell.latvec.e11 * R[iw1][iw2][i].x
						+ GlobalC::ucell.latvec.e21 * R[iw1][iw2][i].y + GlobalC::ucell.latvec.e31 * R[iw1][iw2][i].z;
						Rcar[iw1][iw2][i].y = GlobalC::ucell.latvec.e12 * R[iw1][iw2][i].x
						+ GlobalC::ucell.latvec.e22 * R[iw1][iw2][i].y + GlobalC::ucell.latvec.e32 * R[iw1][iw2][i].z;
						Rcar[iw1][iw2][i].z = GlobalC::ucell.latvec.e13 * R[iw1][iw2][i].x
						+ GlobalC::ucell.latvec.e23 * R[iw1][iw2][i].y + GlobalC::ucell.latvec.e33 * R[iw1][iw2][i].z;
					}
				}
			}

			overlap = new std::complex<double> ***[GlobalV::NLOCAL];
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				overlap[iw1] = new std::complex<double>**[GlobalV::NLOCAL];
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					overlap[iw1][iw2] = new std::complex<double>*[dim];
					for(int ig=0; ig<dim; ig++)
					{
						overlap[iw1][iw2][ig] = new std::complex<double>[NR];
						ModuleBase::GlobalFunc::ZEROS( overlap[iw1][iw2][ig], NR);
					}
				}
			}

			for(int iq=start_q;iq< (start_q + interval_q * nq); iq=iq+interval_q)
			{
				//--------------------------------
				// readin overlap_q
				//--------------------------------
				//std::cout << "begin overlap"<<std::endl;
				std::stringstream ssq;
				ssq << GlobalV::global_out_dir <<"q_"<<icount;
				std::ifstream ifso(ssq.str().c_str());
				if (!ifso)
				{
					ModuleBase::WARNING_QUIT("chi0_hilbert", "Can't find the overlap_q file!");
				}
				ifso >> tmp; //std::cout << tmp << std::endl;
				for(int ig=0; ig<dim; ig++)
				{
					for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
					{
						for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
						{
							for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
							{
								ifso >> overlap[iw1][iw2][ig][ir];
								//std::cout<<"overlap["<<iw1<<"]["<<iw2<<"]["<<ig<<"]["<<ir<<"] = "<<overlap[iw1][iw2][ig][ir] <<std::endl;
							}
						}
					}
				}
				ifso.close();

				time_t start,end;
				start = time(NULL);
				Cal_Chi0s(iq);  // calculate chi0s
				end = time(NULL);
				std::cout<<"chi0s time:: "<<end-start<<" s"<<std::endl;

				time_t start1,end1;
				start1 = time(NULL);
				Cal_T();
				end1 = time(NULL);
				std::cout<<"T time:: "<<end1-start1<<" s"<<std::endl;

				time_t start2,end2;
				start2 = time(NULL);
				Cal_Chi0();   // calculate  chi0
				end2 = time(NULL);
				std::cout<<"chi0 time:: "<<end2-start2<<" s"<<std::endl;

				if(out_chi0)
				{
					plot_chi0(iq);
				}
				for(int i=0; i<nomega; i++)
				{
					std::cout <<"iq = "<<iq<<" chi0["<<i<<"] = "<<chi0[0][i].real()<<"  "<<chi0[0][i].imag()<<std::endl;
				}

				Cal_kernel(iq);

				time_t start3,end3;
				start3 = time(NULL);
				Cal_Chi(iq);  // calculate chi
				end3 = time(NULL);
				std::cout<<"chi time:: "<<end3-start3<<" s"<<std::endl;

				icount++;
			}

			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					delete[] R[iw1][iw2];
				}
				delete[] R[iw1];
			}
			delete[] R;

			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					delete[] Rcar[iw1][iw2];
				}
				delete[] Rcar[iw1];
			}
			delete[] Rcar;

			for(int iw=0; iw<GlobalV::NLOCAL; iw++)
			{
				delete[] Rmax[iw];
			}
			delete[] Rmax;

			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
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
			ModuleBase::WARNING_QUIT("chi0_hilbert","pw, lcao_in_pw, or lcao must be chosen!!!!!!");
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
			std::cout<<"chi0s time:: "<<end-start<<" s"<<std::endl;

			time_t start1, end1;
			start1 = time(NULL);
			Cal_T();
			end1 = time(NULL);
			std::cout<<"T time:: "<<end1-start1<<" s"<<std::endl;

			time_t start2, end2;
			start2 = time(NULL);
			Cal_Chi0();
			end2 = time(NULL);
			std::cout<<"chi time:: "<<end2-start2<<" s"<<std::endl;

			if(out_chi0)
			{
				plot_chi0(iq);
			}
			for(int i=0; i<nomega; i++)
			{
				std::cout <<"iq = "<<iq<<" chi0["<<i<<"] = "<<chi0[0][i].real()<<"  "<<chi0[0][i].imag()<<std::endl;
			}

			Cal_kernel_2D(iq);

			time_t start3, end3;
			start3 = time(NULL);
			Cal_Chi_surface(iq);
			end3 = time(NULL);
			std::cout<<"chi time:: "<<end3-start3<<" s"<<std::endl;
		}
	}
	else
	{
		ModuleBase::WARNING_QUIT("chi0_hilbert","system must be bulk or surface!!!!!!!");
	}

	Delete();

	return;
}

void Chi0_hilbert::Init()
{
	ModuleBase::TITLE("Chi0_hilbert","Init");
	//std::cout << "nbands(init) = " <<GlobalV::NBANDS <<std::endl;
	//std::cout << "oband = " <<oband <<std::endl;
	//std::cout << "nrxx = "<<GlobalC::pw.nrxx<<std::endl;
	//std::cout << "ngmc = " << GlobalC::pw.ngmc <<std::endl;
	//std::cout << "GlobalC::kv.ngk[0] = "<<GlobalC::kv.ngk[0]<<std::endl;
	//std::cout << "dim = "<<dim<<std::endl;
	//std::cout << "GlobalC::kv.nks = " <<GlobalC::kv.nks <<std::endl;
	//std::cout << "GlobalC::kv.nmp = " << GlobalC::kv.nmp[0]<<"  "<<GlobalC::kv.nmp[1]<<"  "<<GlobalC::kv.nmp[2]<<std::endl;

	b_core = new std::complex<double>[GlobalC::pw.ngmc];

	b_summary = new std::complex<double>[GlobalC::pw.ngmc_g];

	b_order = new std::complex<double>[GlobalC::pw.ngmc_g];

	/*psi_r = new std::complex<double>**[GlobalC::kv.nks];
	for(int iq=0; iq<GlobalC::kv.nks; iq++)
	{
		psi_r[iq] = new std::complex<double>*[GlobalV::NBANDS];
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			psi_r[iq][ib] = new std::complex<double>[GlobalC::pw.nrxx];
		}
	}*/

	if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
	{
		psi_r1 = new std::complex<double>*[GlobalV::NBANDS];
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			psi_r1[ib] = new std::complex<double>[GlobalC::pw.nrxx];
		}

		psi_r2 = new std::complex<double>*[GlobalV::NBANDS];
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			psi_r2[ib] = new std::complex<double>[GlobalC::pw.nrxx];
		}
	}
	//std::cout << "psi1 OK" <<std::endl;
	//std::cout << "psi2 OK" <<std::endl;
	std::cout << "psi OK" <<std::endl;

	b = new std::complex<double>**[dim];
	for(int g0=0; g0<dim; g0++)
	{
		b[g0] = new std::complex<double>*[oband];
		for(int ib=0; ib<oband; ib++)
		{
			b[g0][ib] = new std::complex<double>[GlobalV::NBANDS];
		}
	}
	std::cout << "b OK" <<std::endl;

	chi0s = new std::complex<double>*[dim*dim];       // calculate the chi0s (imagniary part of chi0)
	for(int g0=0; g0<dim*dim; g0++)
	{
		chi0s[g0] = new std::complex<double>[nomega];
	}
	std::cout << "chi0s OK" <<std::endl;

	chi0 = new std::complex<double>*[dim*dim];
	for(int g0=0; g0<dim*dim; g0++)
	{
		chi0[g0] = new std::complex<double>[nomega];
	}
	std::cout << "chi0 OK" <<std::endl;

	chi0_gg = new std::complex<double>*[dim];       // split the nomega of chi0
	for(int g0=0; g0<dim; g0++)
	{
		chi0_gg[g0] = new std::complex<double> [dim];
	}
	std::cout << "chi0_gg OK" <<std::endl;

	T = new std::complex<double>*[nomega];              // tranaslation  matrix
	for(int i=0; i<nomega; i++)
	{
		T[i] = new std::complex<double>[nomega];
	}
	std::cout << "T OK" <<std::endl;

	rpa = new std::complex<double>*[dim];
	for(int g0=0; g0<dim; g0++)
	{
		rpa[g0] = new std::complex<double> [dim];
	}
	std::cout << "rpa OK" <<std::endl;

	std::cout <<"system = "<<system<<std::endl;
	if(system == "surface")
	{
		chi0_para = new std::complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			chi0_para[g0] = new std::complex<double> [dim_para];
		}
		std::cout << "chi0_para OK" <<std::endl;

		kernel = new std::complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			kernel[g0] = new std::complex<double> [dim_para];
		}
		std::cout << "kernel OK" <<std::endl;

		rpa1 = new std::complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			rpa1[g0] = new std::complex<double> [dim_para];
		}
		std::cout << "rpa1 OK" <<std::endl;

		chi_para = new std::complex<double>*[dim_para];
		for(int g0=0; g0<dim_para; g0++)
		{
			chi_para[g0] = new std::complex<double> [dim_para];
		}
		std::cout << "chi_para" <<std::endl;
	}
	else
	{
		kernel = new std::complex<double>*[dim];
		for(int g0=0; g0<dim; g0++)
		{
			kernel[g0] = new std::complex<double> [dim];
		}
		std::cout << "kernel OK" <<std::endl;

		chi = new std::complex<double>*[dim];
		for(int g0=0; g0<dim; g0++)
		{
			chi[g0] = new std::complex<double> [dim];
		}
		std::cout << "chi OK" <<std::endl;
	}
	std::cout << "ALL OK"<<std::endl;

	return;
}

void Chi0_hilbert::Delete()
{
	ModuleBase::TITLE("Chi0_hilbert","Delete");
	if(init_finish)
	{
		delete[] b_core;
		delete[] b_summary;
		delete[] b_order;

		for(int i=0; i<GlobalC::kv.nks; i++)
		{
			delete[] cweight[i];
		}
		delete[] cweight;

		if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
		{
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

		std::cout <<"delete OK"<<std::endl;
	}

	return;
}

//----------------------------------------------------------
//    parallel in G     collection and sorting
//----------------------------------------------------------
void Chi0_hilbert::Parallel_G()
{
	ModuleBase::TITLE("Chi0_hilbert","Parallel_G");
	//----------------------------
	// init
	//----------------------------
	num_G_core = new int[GlobalV::DSIZE];
	num_G_dis = new int[GlobalV::DSIZE];
	G_r_core = new double[GlobalC::pw.ngmc];
	num_Gvector_core = new int[GlobalV::DSIZE];
	num_Gvector_dis = new int[GlobalV::DSIZE];
	G_r = new double[GlobalC::pw.ngmc_g];
	Gvec_core = new double[3*GlobalC::pw.ngmc];
	Gvec = new double[3*GlobalC::pw.ngmc_g];
	all_gcar = new ModuleBase::Vector3<double>[GlobalC::pw.ngmc_g];
	flag = new int[GlobalC::pw.ngmc_g];

	for(int i=0;i<GlobalC::pw.ngmc_g;i++)
	{
		flag[i] = i;
	}
	
	ModuleBase::GlobalFunc::ZEROS( num_G_dis, GlobalV::DSIZE);
	ModuleBase::GlobalFunc::ZEROS( num_G_core, GlobalV::DSIZE);
	ModuleBase::GlobalFunc::ZEROS( num_Gvector_dis, GlobalV::DSIZE);
	ModuleBase::GlobalFunc::ZEROS( num_Gvector_core, GlobalV::DSIZE);
	
#ifdef __MPI
	MPI_Allgather( &GlobalC::pw.ngmc, 1, MPI_INT, num_G_core, 1, MPI_INT, POOL_WORLD);
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

	for(int g0=0;g0<GlobalC::pw.ngmc; g0++)
	{
		G_r_core[g0] = GlobalC::pw.get_NormG_cartesian(g0);
		Gvec_core[3*g0] = GlobalC::pw.get_G_cartesian_projection(g0 , 0);
		Gvec_core[3 * g0 + 1] = GlobalC::pw.get_G_cartesian_projection(g0, 1);
		Gvec_core[3 * g0 + 2] = GlobalC::pw.get_G_cartesian_projection(g0, 2);
	}

#ifdef __MPI
	MPI_Allgatherv( G_r_core, GlobalC::pw.ngmc, MPI_DOUBLE, G_r, num_G_core, num_G_dis, MPI_DOUBLE, POOL_WORLD);
	MPI_Allgatherv( Gvec_core, 3*GlobalC::pw.ngmc, MPI_DOUBLE, Gvec, num_Gvector_core, num_Gvector_dis, MPI_DOUBLE, POOL_WORLD);
#endif

	double t1; int t2;
	for(int i=0;i<GlobalC::pw.ngmc_g;i++)
	{
		for(int j=0;j<GlobalC::pw.ngmc_g-i-1;j++)
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

	for(int i=0;i<GlobalC::pw.ngmc_g;i++)
	{
		all_gcar[i].x = Gvec[3*i]; all_gcar[i].y = Gvec[3*i+1]; all_gcar[i].z = Gvec[3*i+2];
		//std::cout<<"all_gcar["<<i<<"]= "<<all_gcar[i].x<<" "<<all_gcar[i].y<<" "<<all_gcar[i].z<<std::endl;
	}

	return;
}

void Chi0_hilbert::Cal_Psi(int iq, std::complex<double> **psi_r)
{
	double phase_x, phase_xy, phase_xyz;
	std::complex<double> exp_tmp;
	for(int ib = 0; ib < GlobalV::NBANDS; ib++)
	{
		ModuleBase::GlobalFunc::ZEROS( GlobalC::UFFT.porter, (GlobalC::pw.nrxx) );
		for(int ig = 0; ig < GlobalC::kv.ngk[iq] ; ig++)
		{
			GlobalC::UFFT.porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(iq,ig)] ] = GlobalC::wf.evc[iq](ib,ig);
		}

		GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter,1);
		int ir=0;
		for(int ix=0; ix<GlobalC::pw.ncx; ix++)
		{
			phase_x = GlobalC::kv.kvec_d[iq].x*ix/GlobalC::pw.ncx;
			for(int iy=0; iy<GlobalC::pw.ncy; iy++)
			{
				phase_xy = phase_x + GlobalC::kv.kvec_d[iq].y*iy/GlobalC::pw.ncy;
				for(int iz=GlobalC::pw.nczp_start; iz<GlobalC::pw.nczp_start+GlobalC::pw.nczp; iz++)
				{
					phase_xyz = (phase_xy + GlobalC::kv.kvec_d[iq].z*iz/GlobalC::pw.ncz) *ModuleBase::TWO_PI;
					exp_tmp = std::complex<double>( cos(phase_xyz), sin(phase_xyz) );
					psi_r[ib][ir] = GlobalC::UFFT.porter[ir]*exp_tmp;
					ir++;
				}

			}
		}
	}

	return;
}


void Chi0_hilbert::Cal_Psi_down(int iq, std::complex<double> **psi_r)
{
	double phase_x, phase_xy, phase_xyz;
	std::complex<double> exp_tmp;
	for(int ib = 0; ib < GlobalV::NBANDS; ib++)
	{
		ModuleBase::GlobalFunc::ZEROS( GlobalC::UFFT.porter, (GlobalC::pw.nrxx) );
		for(int ig = GlobalC::wf.npwx; ig < GlobalC::wf.npwx + GlobalC::kv.ngk[iq] ; ig++)
		{
			GlobalC::UFFT.porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(iq,ig - GlobalC::wf.npwx)] ] = GlobalC::wf.evc[iq](ib,ig);
		}

		GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter,1);
		int ir=0;
		for(int ix=0; ix<GlobalC::pw.ncx; ix++)
		{
			phase_x = GlobalC::kv.kvec_d[iq].x*ix/GlobalC::pw.ncx;
			for(int iy=0; iy<GlobalC::pw.ncy; iy++)
			{
				phase_xy = phase_x + GlobalC::kv.kvec_d[iq].y*iy/GlobalC::pw.ncy;
				for(int iz=GlobalC::pw.nczp_start; iz<GlobalC::pw.nczp_start+GlobalC::pw.nczp; iz++)
				{
					phase_xyz = (phase_xy + GlobalC::kv.kvec_d[iq].z*iz/GlobalC::pw.ncz) *ModuleBase::TWO_PI;
					exp_tmp = std::complex<double>( cos(phase_xyz), sin(phase_xyz) );
					psi_r[ib][ir] = GlobalC::UFFT.porter[ir]*exp_tmp;
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
#ifdef __LCAO
void Chi0_hilbert::Cal_lcao_psi()
{
}

//---------------------------------------------------------------------------
//  calculate the b matrix
//---------------------------------------------------------------------------
void Chi0_hilbert::Cal_b_lcao(int iq, int ik, int iqk)
{
	ModuleBase::TITLE("Chi0_hilbert","Cal_b_lcao");
	double arg;
	ModuleBase::Vector3<double> qk;
	qk = GlobalC::kv.kvec_c[iq] + GlobalC::kv.kvec_c[ik];
	//std::cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<std::endl;

	for(int ig=0; ig<dim; ig++)
	{
		for(int ib1=0; ib1<oband; ib1++)
		{
			ModuleBase::GlobalFunc::ZEROS( b[ig][ib1], GlobalV::NBANDS);
		}
	}

	std::vector<std::vector<std::vector<std::complex<double>>>> phase(GlobalV::NLOCAL,
		std::vector<std::vector<std::complex<double>>>(GlobalV::NLOCAL,
			std::vector<std::complex<double>>(NR)));	// Peize Lin change ptr to std::vector at 2020.01.31

	int MR = 0;
	for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
	{
		for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
		{
			for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
			{
				arg = qk * Rcar[iw1][iw2][ir] * ModuleBase::TWO_PI;
				phase[iw1][iw2][ir] = std::complex<double>( cos(arg),  sin(arg) );
				MR++;
				//std::cout << "phase["<<iw1<<"]["<<iw2<<"]["<<ir<<"] = "<<phase[iw1][iw2][ir]<<std::endl;
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
			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				//------------------------------
				// calculate the matrix
				//------------------------------
				for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
					{
						for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
						{
							b[ig][ib1][ib2] += conj(this->wfc_k_grid[ik][ib1][iw1]) * this->wfc_k_grid[iqk][ib2][iw2] * overlap[iw1][iw2][ig][ir] * phase[iw1][iw2][ir];
						}
					}
				}
			}
		}
	}*/

	std::vector<std::vector<std::complex<double>>> left(oband,
		std::vector<std::complex<double>>(MR));				// Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<std::vector<std::complex<double>>> right(MR,
		std::vector<std::complex<double>>(GlobalV::NBANDS));			// Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<std::vector<std::complex<double>>> b_single(oband,
		std::vector<std::complex<double>>(GlobalV::NBANDS)); 			// Peize Lin change ptr to std::vector at 2020.01.31
	for(int ig=0; ig<dim; ig++)
	{
		for(int ib1=0; ib1<oband; ib1++)
		{
			int r=0;
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
					{
						left[ib1][r] = conj(this->wfc_k_grid[ik][ib1][iw1]) * overlap[iw1][iw2][ig][ir] * phase[iw1][iw2][ir];
						r++;
					}
				}
			}
		}

		for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
		{
			int r=0;
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
					{
						right[r][ib2] = this->wfc_k_grid[iqk][ib2][iw2];
						r++;
					}
				}
			}
		}

		std::complex<double> alpha=1;
		std::complex<double> beta=0;
		char transa = 'N';
		char transb = 'N';
		int M = GlobalV::NBANDS;
		int N = oband;
		int K = MR;
		int lda =GlobalV::NBANDS;
		int ldb = MR;
		int ldc =GlobalV::NBANDS;
		zgemm_(&transa, &transb, &M, &N, &K, &alpha, ModuleBase::GlobalFunc::VECTOR_TO_PTR(right[0]), &lda, 
			ModuleBase::GlobalFunc::VECTOR_TO_PTR(left[0]), &ldb, &beta,  ModuleBase::GlobalFunc::VECTOR_TO_PTR(b_single[0]), &ldc);
		
		for(int ib1=0; ib1<oband; ib1++)
		{
			for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
			{
				b[ig][ib1][ib2] = b_single[ib1][ib2];
			}
		}
	}

	return;
}

#endif

void Chi0_hilbert::Cal_b(int iq, int ik, int iqk, int ispin)
{
	ModuleBase::TITLE("Chi0_hilbert","Cal_b");
	ModuleBase::Vector3<double> qk;
	qk = GlobalC::kv.kvec_c[iq] + GlobalC::kv.kvec_c[ik];
	//std::cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<std::endl;
	double phase_x, phase_xy, phase_xyz;
	ModuleBase::Vector3<double> q = GlobalC::kv.kvec_d[iq];
	std::complex<double> exp_tmp;

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
		for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
		{
			int ir=0;
			for(int ix=0; ix<GlobalC::pw.ncx; ix++)
			{
				phase_x = q.x*ix/GlobalC::pw.ncx;
				for(int iy=0; iy<GlobalC::pw.ncy; iy++)
				{
					phase_xy = phase_x + q.y*iy/GlobalC::pw.ncy;
					for(int iz=GlobalC::pw.nczp_start; iz<GlobalC::pw.nczp_start+GlobalC::pw.nczp; iz++)
					{
						phase_xyz = (phase_xy + q.z*iz/GlobalC::pw.ncz) *ModuleBase::TWO_PI;
						exp_tmp = std::complex<double>(cos(-phase_xyz), sin(-phase_xyz));
						GlobalC::UFFT.porter[ir] = conj(psi_r1[ib1][ir]) * psi_r2[ib2][ir] *exp_tmp;
						ir++;
					}
				}
			}

			GlobalC::pw.FFT_chg.FFT3D( GlobalC::UFFT.porter, -1);
			//for(int g0=0; g0<dim; g0++)
			//{
			//	b[g0][ib1][ib2] = GlobalC::UFFT.porter[ GlobalC::pw.ig2fftc[g0] ];
			//}
			for(int g0=0;g0<GlobalC::pw.ngmc; g0++)
			{
				b_core[g0] = GlobalC::UFFT.porter[ GlobalC::pw.ig2fftc[g0] ];
			}

#ifdef __MPI
			MPI_Allgatherv( b_core, GlobalC::pw.ngmc, mpicomplex, b_summary, num_G_core, num_G_dis, mpicomplex, POOL_WORLD);
#endif

			for(int i=0;i<GlobalC::pw.ngmc_g;i++)
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
	//		for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
	//		{
	//			b[i][ib1][ib2] = std::complex<double>(0.0,0.0);
	//		}
	//	}

	//	for(int ib1=0; ib1<oband; ib1++)
	//	{
	//		for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
	//		{
	//			for(int ig =0; ig< GlobalC::kv.ngk[ik]; ig++)
	//			{
	//				b[ig][ib1][ib2] += conj(GlobalC::wf.evc[ik](ib1,ig)) * sqrt(((GlobalC::kv.kvec_c[ik]+GlobalC::pw.gcar[ig])*(-ModuleBase::TWO_PI/GlobalC::ucell.lat0)).norm2()) * GlobalC::wf.evc[ik](ib2,ig);
	//			}
	//			b[ig][ib1][ib2] /=(GlobalC::wf.ekb[ik][ib1] - GlobalC::wf.ekb[ik][ib2]);
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

        ModuleBase::Vector3<double> qk;
        qk = GlobalC::kv.kvec_c[iq] + GlobalC::kv.kvec_c[ik];
        //std::cout <<"qk = "<<qk.x<<" "<<qk.y<<" "<<qk.z<<std::endl;

        double phase_x, phase_xy, phase_xyz;
        ModuleBase::Vector3<double> q = GlobalC::kv.kvec_d[iq];
        std::cout <<"q_in = "<<q.x<<" "<<q.y<<" "<<q.z<<std::endl;
        std::complex<double> exp_tmp;
        Cal_Psi1(ik);
        Cal_Psi2(iqk);

        for(int ib1=0;ib1<oband; ib1++)
            for(int ib2=0;ib2<GlobalV::NBANDS; ib2++)
                for(int ig=0;ig<dim;ig++)
                {
                    b[ig][ib1][ib2].real() = 0.0;
                    b[ig][ib1][ib2].imag() = 0.0;
                }

        //for(int ib1=0; ib1<oband; ib1++)
        //{
        //        for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
        //        {

        //            for(int ig=0; ig<dim; ig++)
        //            {
        //                int ir=0;
        //                for(int ix=0; ix<GlobalC::pw.ncx; ix++)
        //                {
        //                        phase_x = (q.x + GlobalC::pw.gdirect[ig].x)*ix/GlobalC::pw.ncx;
        //                        for(int iy=0; iy<GlobalC::pw.ncy; iy++)
        //                        {
        //                                phase_xy = phase_x + (q.y + GlobalC::pw.gdirect[ig].y)*iy/GlobalC::pw.ncy;
        //                                for(int iz=GlobalC::pw.nczp_start; iz<GlobalC::pw.nczp_start+GlobalC::pw.nczp; ++iz)
        //                                {
        //                                        phase_xyz = (phase_xy + (q.z + GlobalC::pw.gdirect[ig].z)*iz/GlobalC::pw.ncz) *ModuleBase::TWO_PI;
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
            for(int ix=0;ix<GlobalC::pw.ncx;ix++)
            {
                phase_x = (q.x + GlobalC::pw.gdirect[ig].x)*ix/GlobalC::pw.ncx;
                for(int iy=0; iy<GlobalC::pw.ncy; iy++)
                {
                    phase_xy = phase_x + (q.y + GlobalC::pw.gdirect[ig].y)*iy/GlobalC::pw.ncy;
                    for(int iz=GlobalC::pw.nczp_start; iz<GlobalC::pw.nczp_start+GlobalC::pw.nczp; ++iz)
                    {
                        phase_xyz = (phase_xy + (q.z + GlobalC::pw.gdirect[ig].z)*iz/GlobalC::pw.ncz) *ModuleBase::TWO_PI;
                        exp_tmp.real() = cos(-phase_xyz);
                        exp_tmp.imag() = sin(-phase_xyz);
                        for(int ib1=0;ib1<oband;ib1++)
                        {
                            for(int ib2=0;ib2<GlobalV::NBANDS;ib2++)
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
            for(int ib2=0;ib2<GlobalV::NBANDS; ib2++)
                for(int ig=0;ig<dim;ig++)
                {
                    b[ig][ib1][ib2] = b[ig][ib1][ib2]/GlobalC::pw.ncxyz;
                }

        return;

}*/

// calculate the chi0s matrix
void Chi0_hilbert::Cal_Chi0s(int iq)
{
	ModuleBase::TITLE("Chi0_hilbert","Cal_Chi0s");
	double delta_e,e1,e2;
	std::complex<double> weight1, weight2;
	for(int g0=0; g0<dim; g0++)
	{
		for(int g1=0; g1<dim; g1++)
		{
			ModuleBase::GlobalFunc::ZEROS( chi0s[g1+g0*dim], nomega);
		}
	}

	//---------------------------------------------------------------------------------
	// some test(use our program's wavefunction and other program's occupation to check)
	//---------------------------------------------------------------------------------
	/*double Q[GlobalC::kv.nks][GlobalV::NBANDS];
	std::ifstream ifs("band.dat");
	if(!ifs)
	{
		ModuleBase::WARNING_QUIT("chi0_hilbert","Can't find band.dat");
	}
	for(int i=0;i<GlobalC::kv.nks;i++)
	{
		for(int j=0;j<GlobalV::NBANDS;j++)
		{
			ifs >> Q[i][j];
		}
	}*/

	//---------------------------------------------------------------------------------
	// no spin
	//---------------------------------------------------------------------------------
	if(GlobalV::NSPIN == 1)
	{
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			int iqk = Cal_iq(ik, iq, GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
			if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
			{
				Cal_b(iq, ik, iqk, 0);
			}
			else
			{
				Cal_b_lcao(iq, ik, iqk);
			}

			for(int ib1=0; ib1<oband; ib1++)
			{
				for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
				{
					delta_e = GlobalC::wf.ekb[iqk][ib2] - GlobalC::wf.ekb[ik][ib1];
					//delta_e = Q[ik][ib2] - Q[ik][ib1];
					if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						weight2 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/GlobalC::ucell.omega, 0.0);
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
						weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
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
	else if(GlobalV::NSPIN == 2)
	{
		// spin up
		for(int ik=0; ik<GlobalC::kv.nks/2; ik++)
		{
			int iqk = Cal_iq(ik, iq, GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
			if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
			{
				Cal_b(iq, ik, iqk, 0);
			}
#ifdef __LCAO
			else
			{
				Cal_b_lcao(iq, ik, iqk);
			}
#endif

			for(int ib1=0; ib1<oband; ib1++)
			{
				for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
				{
					delta_e = GlobalC::wf.ekb[iqk][ib2] - GlobalC::wf.ekb[ik][ib1];
					if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						//weight1 = std::complex<double>( (GlobalC::wf.wg(ik,ib1) - GlobalC::wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						//weight2 = std::complex<double>( (GlobalC::wf.wg(ik,ib1) - GlobalC::wf.wg(iqk,ib2)) * (delta_e - e1)/domega/GlobalC::ucell.omega, 0.0);
						weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						weight2 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/GlobalC::ucell.omega, 0.0);
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
						//weight1 = std::complex<double>( (GlobalC::wf.wg(ik,ib1) - GlobalC::wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
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
		for(int ik=GlobalC::kv.nks/2 + 1; ik<GlobalC::kv.nks; ik++)
		{
			int iqk = Cal_iq(ik-GlobalC::kv.nks/2, iq, GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]) + GlobalC::kv.nks/2;
			if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
			{
				Cal_b(iq+GlobalC::kv.nks/2, ik, iqk, 0);
			}
			else
			{
				Cal_b_lcao(iq+GlobalC::kv.nks/2, ik, iqk);
			}
			for(int ib1=0; ib1<oband; ib1++)
			{
				for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
				{
					delta_e = GlobalC::wf.ekb[iqk][ib2] - GlobalC::wf.ekb[ik][ib1];
					if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
					{
						int n = int(delta_e/domega);
						e1 = double(n) * domega;
						e2 = double(n+1) * domega;
						//weight1 = std::complex<double>( (GlobalC::wf.wg(ik,ib1) - GlobalC::wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						//weight2 = std::complex<double>( (GlobalC::wf.wg(ik,ib1) - GlobalC::wf.wg(iqk,ib2)) * (delta_e - e1)/domega/GlobalC::ucell.omega, 0.0);
						weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						weight2 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/GlobalC::ucell.omega, 0.0);
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
						//weight1 = std::complex<double>( (GlobalC::wf.wg(ik,ib1) - GlobalC::wf.wg(iqk,ib2)) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
						weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
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
	else if(GlobalV::NSPIN == 4)
	{
		std::cout<<"NSPIN = "<<GlobalV::NSPIN<<std::endl;
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			int iqk = Cal_iq(ik, iq, GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
			for(int ispin =0; ispin<2; ispin++)
			{
				if(GlobalV::BASIS_TYPE == "pw" || GlobalV::BASIS_TYPE == "lcao_in_pw")
				{
					Cal_b(iq, ik, iqk, ispin);
				}
				else
				{
					Cal_b_lcao(iq, ik, iqk);
				}
				std::cout<<"ik = "<<ik<<" ispin = "<<ispin<<std::endl;
				for(int g=0; g<dim; g++)
				{
					std::cout<<"b["<<g<<"][2][3]"<<" "<<b[g][2][3]<<std::endl;
				}

				for(int ib1=0; ib1<oband; ib1++)
				{
					for(int ib2=0; ib2<GlobalV::NBANDS; ib2++)
					{
						delta_e = GlobalC::wf.ekb[iqk][ib2] - GlobalC::wf.ekb[ik][ib1];
						//delta_e = Q[ik][ib2] - Q[ik][ib1];
						if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
						{
							int n = int(delta_e/domega);
							e1 = double(n) * domega;
							e2 = double(n+1) * domega;
							weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
							weight2 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (delta_e - e1)/domega/GlobalC::ucell.omega, 0.0);
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
							weight1 = std::complex<double>( (cweight[ik][ib1] - cweight[iqk][ib2]) * (e2 - delta_e)/domega/GlobalC::ucell.omega, 0.0);
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
	ModuleBase::TITLE("Chi0_hilbert","Cal_T");
	std::complex<double> M1, M2;
	for(int n1=0; n1<nomega; n1++)
	{
		for(int n=0; n<nomega; n++)
		{
			double n_e = double(n) * domega;
			double n1_e = double(n1) * domega;
			M1 = std::complex<double>( (n_e - n1_e), eta);
			M2 = std::complex<double>( (n_e + n1_e), eta);
			T[n1][n] = 1.0/M1 - 1.0/M2;
		}
	}

	return;
}

// calculate the chi0
void Chi0_hilbert::Cal_Chi0()
{
	ModuleBase::TITLE("Chi0_hilbert","Cal_Chi0");
	std::vector<std::vector<std::complex<double>>> chi0s2(dim*dim,
		std::vector<std::complex<double>>(nomega));			// Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<std::vector<std::complex<double>>> T2(nomega,
		std::vector<std::complex<double>>(nomega));			// Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<std::vector<std::complex<double>>> chi02(dim*dim,
		std::vector<std::complex<double>>(nomega));			// Peize Lin change ptr to std::vector at 2020.01.31
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
	std::complex<double> alpha=1;
	std::complex<double> beta=0;
	char transa = 'N';
	char transb = 'N';
	int M = nomega;
	int N = dim*dim;
	int K = nomega;
	int lda =nomega;
	int ldb =nomega;
	int ldc =nomega;
	
	zgemm_(&transa, &transb, &M, &N, &K, &alpha, ModuleBase::GlobalFunc::VECTOR_TO_PTR(T2[0]), &lda, ModuleBase::GlobalFunc::VECTOR_TO_PTR(chi0s2[0]), &ldb, &beta, ModuleBase::GlobalFunc::VECTOR_TO_PTR(chi02[0]), &ldc);
	
	for(int i=0;i<dim*dim;i++)
	{
		for(int j=0;j<nomega;j++)
		{
			chi0[i][j] = chi02[i][j];
		}
	}

	/*for(int i=0;i<dim*dim;i++)
	{
		ModuleBase::GlobalFunc::ZEROS(chi0[i], nomega);
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

std::complex<double> Chi0_hilbert::f(int k, int j)
{
	std::complex<double> m,n,l;
	m = std::complex<double>( (j - k) * domega, eta);
	n = std::complex<double>( (j + k) * domega, eta);
	l = 1.0/m - 1.0/n;
	return l;
}


// calculate the chi
void Chi0_hilbert::Cal_Chi(int iq)
{
	ModuleBase::TITLE("Chi0_hilbert","Cal_Chi");
	double q = sqrt(((GlobalC::kv.kvec_c[iq])*(ModuleBase::TWO_PI/GlobalC::ucell.lat0)).norm2());
	std::vector<std::complex<double>> gather_chi(nomega);			// Peize Lin change ptr to std::vector at 2020.01.31
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
		//std::cout <<"rpa ok"<<std::endl;

		for(int g0=0; g0<dim; g0++)
		{
			ModuleBase::GlobalFunc::ZEROS( chi[g0], dim);
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
				chi[g0][g1] = std::complex<double>( rpa[g0][g1].real(), rpa[g0][g1].imag());
			}
		}
		std::cout <<"iq = "<<iq<< " epsilon["<<i<<"] = "<<chi[0][0].real()<<"  "<<chi[0][0].imag() <<std::endl;*/

		CMatrixMul(dim, dim, dim, chi0_gg, rpa, chi);
		//std::cout <<"chi ok"<<std::endl;

		gather_chi[i] = chi[0][0];

		std::cout <<"iq = "<<iq<< " chi["<<i<<"] = "<<chi[0][0].real()<<"  "<<chi[0][0].imag() <<std::endl;
		//std::cout << "coulomb_cutoff = " << coulomb_cutoff << std::endl;
		//---------------------------------------------------------------------
		// -Im{epsilon^-1} = -4pi/|q|^2 * Im{chi}
		//---------------------------------------------------------------------
		if(!coulomb_cutoff)
		{
			std::cout <<"iq = "<<iq<<" epsilon["<<i<<"] = "<<8*ModuleBase::PI/q/q*chi[0][0].real()<<"  "<<8*ModuleBase::PI/q/q*chi[0][0].imag()<<std::endl;
		}
		else
		{
			//---------------------------------------------------------------------
			// -Im{epsilon^-1} = -4pi/|q|^2 * (1 - exp^{-q*lz/2}) * Im{chi}
			//---------------------------------------------------------------------
			double factor;
			factor = 1.0 - exp(-q * (GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0)/2.0);
			std::cout <<"iq = "<<iq<<" epsilon["<<i<<"] = "<<8*ModuleBase::PI/q/q * factor * chi[0][0].real()<<"  "<<8*ModuleBase::PI/q/q * factor * chi[0][0].imag()<<std::endl; 
		}

		//--------------------------------------------------------------------
		//  Re{epsilon}  NJP(2014)
		//--------------------------------------------------------------------
		//std::complex<double> M;
		//M = 1.0/(8*ModuleBase::PI/q/q*chi[0][0]);
		//std::cout <<"iq = "<<iq<<" epsMreal["<<i<<"] = "<<M.real()<<std::endl;
		//std::cout <<"iq = "<<iq<<" epsMimag["<<i<<"] = "<<M.imag()<<std::endl;
	}

	if(out_chi)
	{
		std::stringstream sschi;
		sschi << GlobalV::global_out_dir << "chi_"<<iq<<".dat";
		std::ofstream ofschi(sschi.str().c_str());
		for(int i=0;i<nomega;i++)
		{
			ofschi << i * domega <<"   "<<gather_chi[i].real()<<"   "<<gather_chi[i].imag() <<std::endl;
		}
		ofschi.close();
	}

	if(out_epsilon)
	{
		std::stringstream sseps;
		sseps << GlobalV::global_out_dir << "Imeps^-1_"<<iq<<".dat";
		std::ofstream ofseps(sseps.str().c_str());
		ofseps<<"Energy(Ry)"<<"   "<<"-Im{epsilon^-1}"<<std::endl;
		if(!coulomb_cutoff)
		{
			for(int i=0;i<nomega;i++)
			{
				ofseps << i * domega <<"   "<<-8*ModuleBase::PI/q/q*gather_chi[i].imag()<<std::endl;
			}
		}
		else
		{
			double factor;
			factor = 1.0 - exp(-q * (GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0)/2.0);
			for(int i=0;i<nomega;i++)
			{
				ofseps << i * domega <<"   "<<-8*ModuleBase::PI/q/q*factor*gather_chi[i].imag()<<std::endl;
			}			
		}
		ofseps.close();
	}

	return;
}

void Chi0_hilbert::Cal_Chi_surface(int iq)
{
	ModuleBase::TITLE("Chi0_hilbert","Cal_Chi_surface");
	std::complex<double> g;
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
		std::cout <<"iq = "<<iq<< " chi0_para["<<i<<"] = "<<chi0_para[0][0].real()<<"  "<<chi0_para[0][0].imag() <<std::endl;
		Cal_Rpa1();
		std::cout <<"iq = "<<iq<< " epsilon["<<i<<"] = "<<rpa1[0][0].real()<<"  "<<rpa1[0][0].imag() <<std::endl;

		//g = Cal_g(iq);
		//std::cout<<"iq = "<<iq<<" g["<<i<<"] = "<<g.real()<<"  "<<g.imag()<<std::endl;
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
						kernel[g0][g1] = 8.0 * ModuleBase::PI/qg2(iq,g0);
					}
					else
					{
						kernel[g0][g1] = std::complex<double>(0.0, 0.0);
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
				int nz = int(fabs(all_gcar[g0].z * GlobalC::ucell.latvec.e33) + 0.5 );
				sign = (nz%2==0)?1.0:-1.0;
				for(int g1=0; g1<dim; g1++)
				{
					if(g0 == g1)
					{
						kernel[g0][g1] = 8.0 * ModuleBase::PI/qg2(iq,g0) * (1 - sign * exp(-qG(iq,g0)*(GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0)/2.0));
					}
					else
					{
						kernel[g0][g1] = std::complex<double>(0.0, 0.0);
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
					rpa[g0][g1] = -8.0 * ModuleBase::PI/qg2(iq,g0) * chi0_gg[g0][g1];
				}
				else
				{
					rpa[g0][g1] = 1.0 - 8.0 * ModuleBase::PI/qg2(iq,g0) * chi0_gg[g0][g1];
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
			int nz = int(fabs(all_gcar[g0].z * GlobalC::ucell.latvec.e33) + 0.5 );
			sign = (nz%2==0)?1.0:-1.0;
			//std::cout <<"g0= "<<g0<<" nz= "<<nz<<"  sign= "<<sign<<std::endl;
			for(int g1=0; g1<dim; g1++)
			{
				if(g0 != g1)
				{
					rpa[g0][g1] = -8.0 * ModuleBase::PI/qg2(iq,g0) * (1 - sign * exp(-qG(iq,g0)*(GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0)/2.0)) * chi0_gg[g0][g1];
				}
				else
				{
					rpa[g0][g1] = 1.0 - 8.0 * ModuleBase::PI/qg2(iq,g0) * (1 - sign * exp(-qG(iq,g0)*(GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0)/2.0)) * chi0_gg[g0][g1];
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
	/*std::complex<double> a[dim*dim];
	for(int g0=0;g0<dim;g0++)
	{
		for(int g1=0;g1<dim;g1++)
		{
			for(int g1=0;g1<dim;g1++)
		}
	}

	for(int g0=0;g0<dim*dim;g0++)
	{
		std::cout << a[g0]<<"  ";
	}
	std::cout <<std::endl;

	char jobvl = 'V';
	char jobvr = 'V';
	int N = dim;
	int lda = dim;
	std::complex<double>* w = (std::complex<double>*)malloc( sizeof(std::complex<double>) * N);
	int ldvr = dim;
	std::complex<double>* vr = (std::complex<double>*)malloc( sizeof(std::complex<double>) * N * ldvr);
	int ldvl = dim;
	std::complex<double>* vl = (std::complex<double>*)malloc( sizeof(std::complex<double>) * N * ldvl);
	int lwork = N * 4;
	double *rwork = (double*)malloc( sizeof(double) * 2 * N);
	std::complex<double> *work = (std::complex<double>*)malloc( sizeof(std::complex<double>) * lwork);
	int info;

	zgeev_(&jobvl, &jobvr, &N, a, &lda, w, vl , &ldvl , vr, &ldvr, work, &lwork, rwork, &info);

	std::cout <<"eigen  ";
	for(int g0=0;g0<dim;g0++)
	{
		std::cout<<w[g0].real()<<"   "<<w[g0].imag()<<"   ";
	}
	std::cout<<std::endl;*/

	//-------------------------------------------------------------
	// inverse the (I -vx)
	//-------------------------------------------------------------
	int l = Cinv(dim, rpa);
	if (l == 0)
	{
		ModuleBase::WARNING_QUIT("chi0_hilbert","(I-v*chi0) is a singular matrix !!");
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
	qg2 = ((GlobalC::kv.kvec_c[iq]+all_gcar[g0])*(ModuleBase::TWO_PI/GlobalC::ucell.lat0)).norm2();
	
	return qg2;
}

//-------------------------------------------------------------
// inverse a std::complex matrix by Gauss elimination method
//-------------------------------------------------------------
int Chi0_hilbert::Cinv(int n, std::complex<double>** a)
{
	int i,j,k;	//l,u,v,w;
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

		//l=k*n+k;
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

//-------------------------------------------------------------
// Matrix multiplication
//-------------------------------------------------------------
void Chi0_hilbert::CMatrixMul(int m, int n, int l, std::complex<double>** A, std::complex<double>** B, std::complex<double>** C)
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
		std::cout <<"G_direct["<<i<<"] = "<<GlobalC::pw.gdirect[i].x<<" "<<GlobalC::pw.gdirect[i].y<<" "<<GlobalC::pw.gdirect[i].z<<std::endl;
		std::cout <<"flag1["<<i<<"] = "<<flag1[i]<<std::endl;
		std::cout <<"G_para["<<i<<"] = "<<para_g[i][0]<<"  "<<para_g[i][1]<<std::endl;
	}

	return dim1;
}

void Chi0_hilbert::chi0_para_g()
{
	for(int g0=0; g0<dim_para; g0++)
	{
		ModuleBase::GlobalFunc::ZEROS( chi0_para[g0], dim_para);
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
			ModuleBase::GlobalFunc::ZEROS( kernel[g0], dim_para);
		}

		for(int g0=0;g0<dim_para;g0++)
			for(int g1=0;g1<dim_para;g1++)
			{
				if(g0 == g1)
				{
					kernel[g0][g1] = 4* ModuleBase::PI/qg(iq,g0) * (GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0);
				}
				else
				{
					kernel[g0][g1] = std::complex<double>(0.0,0.0);
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
	qG = sqrt((GlobalC::kv.kvec_c[iq].x + all_gcar[g0].x) * (GlobalC::kv.kvec_c[iq].x + all_gcar[g0].x) + (GlobalC::kv.kvec_c[iq].y + all_gcar[g0].y) * (GlobalC::kv.kvec_c[iq].y + all_gcar[g0].y)) * (ModuleBase::TWO_PI/GlobalC::ucell.lat0);
     	
	return qG;
}

double Chi0_hilbert::qg(int iq, int g0)
{
	double qg;
	qg = sqrt((GlobalC::kv.kvec_c[iq].x + para_g[g0][0]) * (GlobalC::kv.kvec_c[iq].x + para_g[g0][0]) + (GlobalC::kv.kvec_c[iq].y + para_g[g0][1]) * (GlobalC::kv.kvec_c[iq].y + para_g[g0][1])) * (ModuleBase::TWO_PI/GlobalC::ucell.lat0);
	
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
		ModuleBase::WARNING_QUIT("chi0_hilbert","(I-v*chi0) is a singular matrix !!");
	}

	return;
}

void Chi0_hilbert::chi_para_g()
{
	CMatrixMul(dim_para, dim_para, dim_para, chi0_para, rpa1, chi_para);

	return;
}


std::complex<double> Chi0_hilbert:: Cal_g(int iq)
{
    std::complex<double> g;
    double L = GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0;
    double dz = L/GlobalC::pw.ncz;
    double q =  sqrt(((GlobalC::kv.kvec_c[iq])*(ModuleBase::TWO_PI/GlobalC::ucell.lat0)).norm2());
    
    g = std::complex<double>(0.0,0.0);

    for(int z0=0;z0<GlobalC::pw.ncz;z0++)
        for(int z1=0;z1<GlobalC::pw.ncz;z1++)
        {
            double exp_phase = exp(-q*(z0+z1)*dz);
            g += chi_para[z0][z1] * exp_phase * dz * dz;
        }
    g = -4*ModuleBase::PI/q * g;

    return g;
}

//----------------------------------------------------
// plot function
//----------------------------------------------------
void Chi0_hilbert::plot_chi0(int iq)
{
	std::stringstream ss;
	ss << GlobalV::global_out_dir <<"chi0"<<"_"<<iq<<".dat";
	std::ofstream ofs(ss.str().c_str());
	for(int i=0;i<nomega;i++)
	{
		ofs <<i*domega<<"    "<<chi0[0][i].real()<<"    "<<chi0[0][i].imag()<<std::endl;
	}
	ofs.close();

	return;
}
