#include "to_wannier90.h"
#include "../src_lcao/global_fp.h" // mohan add 2021-01-30, this module should be modified
 


toWannier90::toWannier90(int num_kpts, Matrix3 recip_lattice)
{
	this->num_kpts = num_kpts;
	this->recip_lattice = recip_lattice;
	if(NSPIN==1 || NSPIN==4) this->cal_num_kpts = this->num_kpts;
	else if(NSPIN==2) this->cal_num_kpts = this->num_kpts/2;

}


toWannier90::~toWannier90()
{
	if(num_exclude_bands > 0) delete[] exclude_bands;
	if(BASIS_TYPE == "lcao") delete[] unk_inLcao;
}


void toWannier90::init_wannier()
{	
	this->read_nnkp();
	
	if(NSPIN == 2)
	{
		wannier_spin = INPUT.wannier_spin;
		if(wannier_spin == "up") start_k_index = 0;
		else if(wannier_spin == "down") start_k_index = num_kpts/2;
		else
		{
			WARNING_QUIT("toWannier90::init_wannier","Error wannier_spin set,is not \"up\" or \"down\" ");
		}
	}
	
	if(BASIS_TYPE == "pw")
	{
		writeUNK(wf.evc);
		outEIG();
		cal_Mmn(wf.evc);
		cal_Amn(wf.evc);
	}
	else if(BASIS_TYPE == "lcao")
	{
		getUnkFromLcao();
		cal_Amn(this->unk_inLcao);
		cal_Mmn(this->unk_inLcao);
		writeUNK(this->unk_inLcao);
		outEIG();
	}

	/*
	if(MY_RANK==0)
	{
		if(BASIS_TYPE == "pw")
		{
			cal_Amn(wf.evc);
			cal_Mmn(wf.evc);
			writeUNK(wf.evc);
			outEIG();
		}
		else if(BASIS_TYPE == "lcao")
		{
			getUnkFromLcao();
			cal_Amn(this->unk_inLcao);
			cal_Mmn(this->unk_inLcao);
			writeUNK(this->unk_inLcao);
			outEIG();
		}
	}
	*/
	
}

void toWannier90::read_nnkp()
{
	// read *.nnkp file
	// 检查 正格矢，倒格矢，k点坐标，试探轨道投影，每个k点的近邻k点，需要排除的能带指标
	
	wannier_file_name = INPUT.NNKP;
	wannier_file_name = wannier_file_name.substr(0,wannier_file_name.length() - 5);

	ofs_running << "reading the " << wannier_file_name << ".nnkp file." << endl;
	
	ifstream nnkp_read(INPUT.NNKP.c_str(), ios::in);
	
	if(!nnkp_read) WARNING_QUIT("toWannier90::read_nnkp","Error during readin parameters.");
	
	if( SCAN_BEGIN(nnkp_read,"real_lattice") )
	{
		Matrix3 real_lattice_nnkp;
		nnkp_read >> real_lattice_nnkp.e11 >> real_lattice_nnkp.e12 >> real_lattice_nnkp.e13
				  >> real_lattice_nnkp.e21 >> real_lattice_nnkp.e22 >> real_lattice_nnkp.e23
				  >> real_lattice_nnkp.e31 >> real_lattice_nnkp.e32 >> real_lattice_nnkp.e33;
				  
		real_lattice_nnkp = real_lattice_nnkp / ucell.lat0_angstrom;
		
		if(abs(real_lattice_nnkp.e11 - ucell.latvec.e11) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e12 - ucell.latvec.e12) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e13 - ucell.latvec.e13) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e21 - ucell.latvec.e21) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e22 - ucell.latvec.e22) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e23 - ucell.latvec.e23) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e31 - ucell.latvec.e31) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e32 - ucell.latvec.e32) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		if(abs(real_lattice_nnkp.e33 - ucell.latvec.e33) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error real_lattice in *.nnkp file");
		
	}
	
	if( SCAN_BEGIN(nnkp_read,"recip_lattice") )
	{
		Matrix3 recip_lattice_nnkp;
		nnkp_read >> recip_lattice_nnkp.e11 >> recip_lattice_nnkp.e12 >> recip_lattice_nnkp.e13
				  >> recip_lattice_nnkp.e21 >> recip_lattice_nnkp.e22 >> recip_lattice_nnkp.e23
				  >> recip_lattice_nnkp.e31 >> recip_lattice_nnkp.e32 >> recip_lattice_nnkp.e33;
		
		const double tpiba_angstrom = TWO_PI / ucell.lat0_angstrom;
		recip_lattice_nnkp = recip_lattice_nnkp / tpiba_angstrom;
		
		if(abs(recip_lattice_nnkp.e11 - ucell.G.e11) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e12 - ucell.G.e12) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e13 - ucell.G.e13) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e21 - ucell.G.e21) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e22 - ucell.G.e22) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e23 - ucell.G.e23) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e31 - ucell.G.e31) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e32 - ucell.G.e32) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
		if(abs(recip_lattice_nnkp.e33 - ucell.G.e33) > 1.0e-4) 
			WARNING_QUIT("toWannier90::read_nnkp","Error recip_lattice in *.nnkp file");
	}
	
	if( SCAN_BEGIN(nnkp_read,"kpoints") )
	{
		int numkpt_nnkp;
		READ_VALUE(nnkp_read, numkpt_nnkp);
		if( (NSPIN == 1 || NSPIN == 4) && numkpt_nnkp != kv.nkstot ) WARNING_QUIT("toWannier90::read_nnkp","Error kpoints in *.nnkp file");
		else if(NSPIN == 2 && numkpt_nnkp != (kv.nkstot/2))	WARNING_QUIT("toWannier90::read_nnkp","Error kpoints in *.nnkp file");
	
		Vector3<double> *kpoints_direct_nnkp = new Vector3<double>[numkpt_nnkp];
		for(int ik = 0; ik < numkpt_nnkp; ik++)
		{
			nnkp_read >> kpoints_direct_nnkp[ik].x >> kpoints_direct_nnkp[ik].y >> kpoints_direct_nnkp[ik].z;
			if(abs(kpoints_direct_nnkp[ik].x - kv.kvec_d[ik].x) > 1.0e-4) 
				WARNING_QUIT("toWannier90::read_nnkp","Error kpoints in *.nnkp file");
			if(abs(kpoints_direct_nnkp[ik].y - kv.kvec_d[ik].y) > 1.0e-4) 
				WARNING_QUIT("toWannier90::read_nnkp","Error kpoints in *.nnkp file");
			if(abs(kpoints_direct_nnkp[ik].z - kv.kvec_d[ik].z) > 1.0e-4) 
				WARNING_QUIT("toWannier90::read_nnkp","Error kpoints in *.nnkp file");
		}
				
		delete[] kpoints_direct_nnkp;
		
		//判断gamma only
		Vector3<double> my_gamma_point(0.0,0.0,0.0);
		//if( (kv.nkstot == 1) && (kv.kvec_d[0] == my_gamma_point) ) gamma_only_wannier = true;
	} 
	
	if(NSPIN!=4)
	{
		if( SCAN_BEGIN(nnkp_read,"projections") )
		{
			READ_VALUE(nnkp_read, num_wannier);
			// test
			//ofs_running << "num_wannier = " << num_wannier << endl;
			// test
			if(num_wannier < 0)
			{
				WARNING_QUIT("toWannier90::read_nnkp","wannier number is lower than 0");
			}
			
			R_centre = new Vector3<double>[num_wannier];
			L = new int[num_wannier];
			m = new int[num_wannier];
			rvalue = new int[num_wannier];
			Vector3<double>* z_axis = new Vector3<double>[num_wannier];
			Vector3<double>* x_axis = new Vector3<double>[num_wannier];
			alfa = new double[num_wannier];
			
			
			for(int count = 0; count < num_wannier; count++)
			{
				nnkp_read >> R_centre[count].x >> R_centre[count].y >> R_centre[count].z;
				nnkp_read >> L[count] >> m[count];
				READ_VALUE(nnkp_read,rvalue[count]);
				nnkp_read >> z_axis[count].x >> z_axis[count].y >> z_axis[count].z;
				nnkp_read >> x_axis[count].x >> x_axis[count].y >> x_axis[count].z;
				READ_VALUE(nnkp_read,alfa[count]);			
			}
			
		}
	}
	else
	{
		WARNING_QUIT("toWannier90::read_nnkp","noncolin spin is not done yet");
	}

	if( SCAN_BEGIN(nnkp_read,"nnkpts") )
	{
		READ_VALUE(nnkp_read, nntot);
		nnlist.resize(kv.nkstot);
		nncell.resize(kv.nkstot);
		for(int ik = 0; ik < kv.nkstot; ik++)
		{
			nnlist[ik].resize(nntot);
			nncell[ik].resize(nntot);
		}
		
		int numkpt_nnkp;
		if(NSPIN == 1 || NSPIN == 4) numkpt_nnkp = kv.nkstot;
		else if(NSPIN == 2) numkpt_nnkp = kv.nkstot/2;
		else throw runtime_error("numkpt_nnkp uninitialized in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		
		for(int ik = 0; ik < numkpt_nnkp; ik++)
		{
			for(int ib = 0; ib < nntot; ib++)
			{
				int ik_nnkp;
				nnkp_read >> ik_nnkp;
				if(ik_nnkp != (ik+1)) WARNING_QUIT("toWannier90::read_nnkp","error nnkpts in *.nnkp file");
				nnkp_read >> nnlist[ik][ib];
				nnkp_read >> nncell[ik][ib].x >> nncell[ik][ib].y >> nncell[ik][ib].z;
				nnlist[ik][ib]--; // this is c++ , begin from 0
			}
			
		}
	}
	
	if( SCAN_BEGIN(nnkp_read,"exclude_bands") )
	{
		READ_VALUE(nnkp_read, num_exclude_bands);
		if(num_exclude_bands > 0) exclude_bands = new int[num_exclude_bands];
		else if(num_exclude_bands < 0) WARNING_QUIT("toWannier90::read_nnkp","the exclude bands is wrong , please check *.nnkp file.");
		
		if(num_exclude_bands > 0)
		{
			for(int i = 0; i < num_exclude_bands; i++)
			{
				READ_VALUE(nnkp_read, exclude_bands[i]);
				exclude_bands[i]--; // this is c++ , begin from 0
			}
		}
	}
	
	// test by jingan
	//ofs_running << "num_exclude_bands = " << num_exclude_bands << endl;
	//for(int i = 0; i < num_exclude_bands; i++)
	//{
	//	ofs_running << "exclude_bands : " << exclude_bands[i] << endl;
	//}
	// test by jingan
	
	nnkp_read.close();
	
	// 设置试探轨道参数
	for(int i = 0; i < num_wannier; i++)
	{
		R_centre[i] = R_centre[i] * ucell.latvec;
		m[i] = m[i] - 1; // ABACUS and wannier90 对磁角动量m的定义不一样，ABACUS是从0开始的，wannier90是从1开始的
	}
	
	// test by jingan
	//ofs_running << "num_wannier is " << num_wannier << endl;
	//for(int i = 0; i < num_wannier; i++)
	//{
	//	ofs_running << "num_wannier" << endl;
	//	ofs_running << L[i] << " " << m[i] << " " << rvalue[i] << " " << alfa[i] << endl;
	//}
	// test by jingan
	
	// 设置exclude_bands
	tag_cal_band = new bool[NBANDS];
	if(NBANDS <= num_exclude_bands) WARNING_QUIT("toWannier90::read_nnkp","you set the band numer is not enough, please add bands number.");
	if(num_exclude_bands == 0)
	{
		for(int ib = 0; ib < NBANDS; ib++) tag_cal_band[ib] = true;
	}
	else
	{
		for(int ib = 0; ib < NBANDS; ib++)
		{
			tag_cal_band[ib] = true;
			for(int ibb = 0; ibb < num_exclude_bands; ibb++)
			{
				if(exclude_bands[ibb] == ib) 
				{
					tag_cal_band[ib] = false;
					break;
				}
			}
		}
	}
	
	if(num_exclude_bands < 0) num_bands = NBANDS;
	else num_bands = NBANDS - num_exclude_bands;
	
	
}

void toWannier90::outEIG()
{
	if(MY_RANK == 0)
	{
		string fileaddress = global_out_dir + wannier_file_name + ".eig";
		ofstream eig_file( fileaddress.c_str() );
		for(int ik = start_k_index; ik < (cal_num_kpts+start_k_index); ik++)
		{
			int index_band = 0;
			for(int ib = 0; ib < NBANDS; ib++)
			{
				if(!tag_cal_band[ib]) continue;
				index_band++;
				eig_file << setw(5) << index_band << setw(5) << ik+1-start_k_index
						 << setw(18) << showpoint << fixed << setprecision(12) 
						 << wf.ekb[ik][ib] * Ry_to_eV << endl;
			}
		}
		
		eig_file.close();
	}
}


void toWannier90::writeUNK(const ComplexMatrix *wfc_pw)
{

	/*
	complex<double> *porter = new complex<double>[pw.nrxx];
	
	for(int ik = start_k_index; ik < (cal_num_kpts+start_k_index); ik++)
	{
		stringstream name;
		if(NSPIN==1 || NSPIN==4)
		{
			name << global_out_dir << "UNK" << setw(5) << setfill('0') << ik+1 << ".1" ;
		}
		else if(NSPIN==2)
		{
			if(wannier_spin=="up") name << global_out_dir << "UNK" << setw(5) << setfill('0') << ik+1-start_k_index << ".1" ;
			else if(wannier_spin=="down") name << global_out_dir << "UNK" << setw(5) << setfill('0') << ik+1-start_k_index << ".2" ;
		}
		
		ofstream unkfile(name.str());
		
		unkfile << setw(12) << pw.ncx << setw(12) << pw.ncy << setw(12) << pw.ncz << setw(12) << ik+1 << setw(12) << num_bands << endl;
		
		for(int ib = 0; ib < NBANDS; ib++)
		{
			if(!tag_cal_band[ib]) continue;
			//complex<double> *porter = UFFT.porter;
			//  u_k in real space
			ZEROS(porter, pw.nrxx);
			for (int ig = 0; ig < kv.ngk[ik]; ig++)
			{
				porter[pw.ig2fftw[wf.igk(ik, ig)]] = wfc_pw[ik](ib, ig);
			}
			pw.FFT_wfc.FFT3D(porter, 1);
			
			for(int k=0; k<pw.ncz; k++)
			{
				for(int j=0; j<pw.ncy; j++)
				{
					for(int i=0; i<pw.ncx; i++)
					{
						if(!gamma_only_wannier)
						{
							unkfile << setw(20) << setprecision(9) << setiosflags(ios::scientific) << porter[i*pw.ncy*pw.ncz + j*pw.ncz + k].real()
									<< setw(20) << setprecision(9) << setiosflags(ios::scientific) << porter[i*pw.ncy*pw.ncz + j*pw.ncz + k].imag() 
									//jingan test
									//<< "       " << setw(12) << setprecision(9) << setiosflags(ios::scientific) << abs(porter[i*pw.ncy*pw.ncz + j*pw.ncz + k])
									<< endl;
						}
						else
						{
							double zero = 0.0;
							unkfile << setw(20) << setprecision(9) << setiosflags(ios::scientific) << abs( porter[i*pw.ncy*pw.ncz + j*pw.ncz + k] )
									<< setw(20) << setprecision(9) << setiosflags(ios::scientific) << zero
									//jingan test
									//<< "       " << setw(12) << setprecision(9) << setiosflags(ios::scientific) << abs(porter[i*pw.ncy*pw.ncz + j*pw.ncz + k])
									<< endl;
						}
					}
				}
			}
			
			
		}

		
		unkfile.close();
		
	}
	
	delete[] porter;
	*/
	
#ifdef __MPI
	// num_z: how many planes on processor 'ip'
	int *num_z = new int[NPROC_IN_POOL];
	ZEROS(num_z, NPROC_IN_POOL);
	for (int iz=0;iz<pw.nbz;iz++)
	{
		int ip = iz % NPROC_IN_POOL;
		num_z[ip] += pw.bz;
	}	

	// start_z: start position of z in 
	// processor ip.
	int *start_z = new int[NPROC_IN_POOL];
	ZEROS(start_z, NPROC_IN_POOL);
	for (int ip=1;ip<NPROC_IN_POOL;ip++)
	{
		start_z[ip] = start_z[ip-1]+num_z[ip-1];
	}	

	// which_ip: found iz belongs to which ip.
	int *which_ip = new int[pw.ncz];
	ZEROS(which_ip, pw.ncz);
	for(int iz=0; iz<pw.ncz; iz++)
	{
		for(int ip=0; ip<NPROC_IN_POOL; ip++)
		{
			if(iz>=start_z[NPROC_IN_POOL-1]) 
			{
				which_ip[iz] = NPROC_IN_POOL-1;
				break;
			}
			else if(iz>=start_z[ip] && iz<start_z[ip+1])
			{
				which_ip[iz] = ip;
				break;
			}
		}
	}
	
	
	// only do in the first pool.
	complex<double> *porter = new complex<double>[pw.nrxx];
	int nxy = pw.ncx * pw.ncy;
	complex<double> *zpiece = new complex<double>[nxy];
	
	if(MY_POOL==0)
	{
		for(int ik = start_k_index; ik < (cal_num_kpts+start_k_index); ik++)
		{
			ofstream unkfile;
			
			if(MY_RANK == 0)
			{
				stringstream name;
				if(NSPIN==1 || NSPIN==4)
				{
					name << global_out_dir << "UNK" << setw(5) << setfill('0') << ik+1 << ".1" ;
				}
				else if(NSPIN==2)
				{
					if(wannier_spin=="up") name << global_out_dir << "UNK" << setw(5) << setfill('0') << ik+1-start_k_index << ".1" ;
					else if(wannier_spin=="down") name << global_out_dir << "UNK" << setw(5) << setfill('0') << ik+1-start_k_index << ".2" ;
				}
				
				unkfile.open(name.str(),ios::out);
				
				unkfile << setw(12) << pw.ncx << setw(12) << pw.ncy << setw(12) << pw.ncz << setw(12) << ik+1 << setw(12) << num_bands << endl;
			}
			
			for(int ib = 0; ib < NBANDS; ib++)
			{
				if(!tag_cal_band[ib]) continue;
				
				ZEROS(porter, pw.nrxx);
				for (int ig = 0; ig < kv.ngk[ik]; ig++)
				{
					porter[pw.ig2fftw[wf.igk(ik, ig)]] = wfc_pw[ik](ib, ig);
				}
				pw.FFT_wfc.FFT3D(porter, 1);

				// save the rho one z by one z.
				for(int iz=0; iz<pw.ncz; iz++)
				{
					// tag must be different for different iz.
					ZEROS(zpiece, nxy);
					int tag = iz;
					MPI_Status ierror;

					// case 1: the first part of rho in processor 0.
					if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
					{
						for(int ir=0; ir<nxy; ir++)
						{
							zpiece[ir] = porter[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
						}
					}
					// case 2: > first part rho: send the rho to 
					// processor 0.
					else if(which_ip[iz] == RANK_IN_POOL )
					{
						for(int ir=0; ir<nxy; ir++)
						{
							zpiece[ir] = porter[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
						}
						MPI_Send(zpiece, nxy, MPI_DOUBLE_COMPLEX, 0, tag, POOL_WORLD);
					}

					// case 2: > first part rho: processor 0 receive the rho
					// from other processors
					else if(RANK_IN_POOL==0)
					{
						MPI_Recv(zpiece, nxy, MPI_DOUBLE_COMPLEX, which_ip[iz], tag, POOL_WORLD, &ierror);
					}

					// write data	
					if(MY_RANK==0)
					{
						for(int iy=0; iy<pw.ncy; iy++)
						{
							for(int ix=0; ix<pw.ncx; ix++)
							{
								unkfile << setw(20) << setprecision(9) << setiosflags(ios::scientific) << zpiece[ix*pw.ncy+iy].real()
										<< setw(20) << setprecision(9) << setiosflags(ios::scientific) << zpiece[ix*pw.ncy+iy].imag() 
										<< endl;
							}
						}
					}
				}// end iz
				MPI_Barrier(POOL_WORLD);
			}
			
			if(MY_RANK == 0)
			{
				unkfile.close();
			}
		
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	delete[] num_z;
	delete[] start_z;
	delete[] which_ip;
	delete[] porter;
	delete[] zpiece;

#endif	
	
}







void toWannier90::cal_Amn(const ComplexMatrix *wfc_pw)
{
	// 第一步：建立实球谐函数lm在某个k点下的平面波基组下的表格（矩阵）	
	// 第二步：将试探轨道的径向部分向某个k点下平面波投影
	// 第三步：获取试探轨道在某个k点下平面波基组下的投影
	const int pwNumberMax = wf.npwx;
	
	ofstream Amn_file;
	
	if(MY_RANK == 0)
	{
		time_t  time_now = time(NULL);
		string fileaddress = global_out_dir + wannier_file_name + ".amn";
		Amn_file.open( fileaddress.c_str() , ios::out);
		Amn_file << " Created on " << ctime(&time_now);
		Amn_file << setw(12) << num_bands << setw(12) << cal_num_kpts << setw(12) << num_wannier << endl;
	}
	
	ComplexMatrix *trial_orbitals = new ComplexMatrix[cal_num_kpts];
	for(int ik = 0; ik < cal_num_kpts; ik++)
	{
		trial_orbitals[ik].create(num_wannier,pwNumberMax);
		produce_trial_in_pw(ik,trial_orbitals[ik]);
	}	
	
	// test by jingan
	//ofs_running << __FILE__ << __LINE__ << "start_k_index = " << start_k_index << "  cal_num_kpts = " << cal_num_kpts << endl;
	// test by jingan

	for(int ik = start_k_index; ik < (cal_num_kpts+start_k_index); ik++)
	{
		for(int iw = 0; iw < num_wannier; iw++)
		{
			int index_band = 0;
			for(int ib = 0; ib < NBANDS; ib++)
			{
				if(!tag_cal_band[ib]) continue;
				index_band++;
				complex<double> amn(0.0,0.0);
				complex<double> amn_tem(0.0,0.0);
				for(int ig = 0; ig < pwNumberMax; ig++)
				{
					int cal_ik = ik - start_k_index;
					amn_tem = amn_tem + conj( wfc_pw[ik](ib,ig) ) * trial_orbitals[cal_ik](iw,ig);
				}
				
				MPI_Allreduce(&amn_tem , &amn , 1, MPI_DOUBLE_COMPLEX , MPI_SUM , POOL_WORLD);
				
				if(MY_RANK == 0)
				{
					Amn_file << setw(5) << index_band << setw(5) << iw+1 << setw(5) << ik+1-start_k_index 
							 << setw(18) << showpoint << fixed << setprecision(12) << amn.real() 
							 << setw(18) << showpoint << fixed << setprecision(12) << amn.imag()
							 //jingan test
							 //<< "   " << setw(18) << setprecision(13) << abs(amn)
							 << endl;
				}
			}
		}
	}
	

	
	if(MY_RANK == 0) Amn_file.close();
	
	delete[] trial_orbitals;
	
}



void toWannier90::cal_Mmn(const ComplexMatrix *wfc_pw)
{	
	// test by jingan
	//ofs_running << __FILE__ << __LINE__ << " cal_num_kpts = " << cal_num_kpts << endl;
	// test by jingan
	
	ofstream mmn_file;
	
	if(MY_RANK == 0)
	{
		string fileaddress = global_out_dir + wannier_file_name + ".mmn";
		mmn_file.open( fileaddress.c_str() , ios::out);	
		
		time_t  time_now = time(NULL);
		mmn_file << " Created on " << ctime(&time_now);
		mmn_file << setw(12) << num_bands << setw(12) << cal_num_kpts << setw(12) << nntot << endl;
	}
	
	/*
	ComplexMatrix Mmn(NBANDS,NBANDS);
	if(gamma_only_wannier)
	{
		for(int ib = 0; ib < nntot; ib++)
		{
			Vector3<double> phase_G = nncell[0][ib];
			for(int m = 0; m < NBANDS; m++)
			{
				if(!tag_cal_band[m]) continue;
				for(int n = 0; n <= m; n++)
				{
					if(!tag_cal_band[n]) continue;
					complex<double> mmn_tem = gamma_only_cal(m,n,wfc_pw,phase_G);
					Mmn(m,n) = mmn_tem;
					if(m!=n) Mmn(n,m) = Mmn(m,n);				
				}
			}
		}
	}
	*/
	
	for(int ik = 0; ik < cal_num_kpts; ik++)
	{
		for(int ib = 0; ib < nntot; ib++)
		{
			int ikb = nnlist[ik][ib];             // ik+b : ik的近邻k点	
			
			Vector3<double> phase_G = nncell[ik][ib];
			
			if(MY_RANK == 0)
			{
				mmn_file << setw(5) << ik+1 << setw(5) << ikb+1 << setw(5) 
						 << int(phase_G.x) << setw(5) << int(phase_G.y) << setw(5) << int(phase_G.z) 
						 << endl;
			}
		
			for(int m = 0; m < NBANDS; m++)
			{
				if(!tag_cal_band[m]) continue;
				for(int n = 0; n < NBANDS; n++)
				{
					if(!tag_cal_band[n]) continue;
					complex<double> mmn(0.0,0.0);
				
					if(!gamma_only_wannier)
					{
						int cal_ik = ik + start_k_index;
						int cal_ikb = ikb + start_k_index;												
						// test by jingan
						//ofs_running << __FILE__ << __LINE__ << "cal_ik = " << cal_ik << "cal_ikb = " << cal_ikb << endl;
						// test by jingan
						//complex<double> *unk_L_r = new complex<double>[pw.nrxx];
						//ToRealSpace(cal_ik,n,wfc_pw,unk_L_r,phase_G);				
						//mmn = unkdotb(unk_L_r,cal_ikb,m,wfc_pw);
						mmn = unkdotkb(cal_ik,cal_ikb,n,m,phase_G,wfc_pw);
						//delete[] unk_L_r;
					}
					else
					{
						//ofs_running << "gamma only test" << endl;
						//mmn = Mmn(n,m);
					}
					
					if(MY_RANK == 0)
					{
						mmn_file << setw(18) << setprecision(12) << showpoint << fixed << mmn.real() 
								 << setw(18) << setprecision(12) << showpoint << fixed << mmn.imag()
								 // jingan test
								 //<< "    " << setw(12) << setprecision(9) << abs(mmn)
								 << endl;				
					}
				}
			}
		}
	
	}
	
	if(MY_RANK == 0) mmn_file.close();
	
}


void toWannier90::produce_trial_in_pw(const int &ik, ComplexMatrix &trial_orbitals_k)
{
	// 检查参数是否正确
	for(int i =0; i < num_wannier; i++)
	{
		if(L[i] < -5 || L[i] > 3) cout << "toWannier90::produce_trial_in_pw() your L angular momentum is wrong , please check !!! " << endl;
	
		if(L[i] >= 0) 
		{
			if(m[i] < 0 || m[i] > 2*L[i]) cout << "toWannier90::produce_trial_in_pw() your m momentum is wrong , please check !!! " << endl;
		}
		else
		{
			if(m[i] < 0 || m[i] > -L[i]) cout << "toWannier90::produce_trial_in_pw() your m momentum is wrong , please check !!! " << endl;
		
		}
	}
	
	const int npw = kv.ngk[ik];
	const int npwx = wf.npwx;
	const int total_lm = 16;
	matrix ylm(total_lm,npw);               //所有类型的球谐函数
	//matrix wannier_ylm(num_wannier,npw);    //要试探轨道的使用的球谐函数
	double bs2, bs3, bs6, bs12;
	bs2 = 1.0/sqrt(2.0);
	bs3 = 1.0/sqrt(3.0);
	bs6 = 1.0/sqrt(6.0);
	bs12 = 1.0/sqrt(12.0);
	
	Vector3<double> *gk = new Vector3<double>[npw];
	for(int ig = 0; ig < npw; ig++)
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);  // k+G矢量
	}
	
	Mathzone::Ylm_Real(total_lm, npw, gk, ylm);
	
	// test by jingan
	//ofs_running << "the mathzone::ylm_real is successful!" << endl;
	//ofs_running << "produce_trial_in_pw: num_wannier is " << num_wannier << endl;
	// test by jingan
	
	
	// 1.生成径向轨道在某个k点平面波基组的投影
	const int mesh_r = 333; 		//描述径向函数所需要的格点数
	const double dx = 0.025; 		//固定间隔，用于生成非固定间隔的dr来提高精度,这个值很巧妙
	const double x_min = -6.0;  	// 用于生成dr和r的起始点
	matrix r(num_wannier,mesh_r);   //不同alfa的径向函数的r
	matrix dr(num_wannier,mesh_r);  //不同alfa的径向函数的每个r点的间隔
	matrix psi(num_wannier,mesh_r); //径向函数psi in 实空间
	matrix psir(num_wannier,mesh_r);// psi * r in 实空间
	matrix psik(num_wannier,npw);   //径向函数在某个k点下倒空间的投影
	
	// 生成r,dr
	for(int i = 0; i < num_wannier; i++)
	{
		double x = 0;
		for(int ir = 0; ir < mesh_r; ir++)
		{
			x = x_min + ir * dx;
			r(i,ir) = exp(x) / alfa[i];
			dr(i,ir) = dx * r(i,ir);
		}
		
	}
	
	// 生成psi
	for(int i = 0; i < num_wannier; i++)
	{
		double alfa32 = pow(alfa[i],3.0/2.0);
		double alfa_new = alfa[i];
		int wannier_index = i;
		
		if(rvalue[i] == 1)
		{
			for(int ir = 0; ir < mesh_r; ir++)
			{
				psi(wannier_index,ir) = 2.0 * alfa32 * exp( -alfa_new * r(wannier_index,ir) );
			}
		}
	
		if(rvalue[i] == 2)
		{
			for(int ir = 0; ir < mesh_r; ir++)
			{
				psi(wannier_index,ir) = 1.0/sqrt(8.0) * alfa32
										* (2.0 - alfa_new * r(wannier_index,ir))
										* exp( -alfa_new * r(wannier_index,ir) * 0.5 );
			}
		}

		if(rvalue[i] == 3)
		{
			for(int ir = 0; ir < mesh_r; ir++)
			{
				psi(wannier_index,ir) = sqrt(4.0/27.0) * alfa32
										* ( 1.0 - 2.0/3.0 * alfa_new * r(wannier_index,ir) + 2.0/27.0 * pow(alfa_new,2.0) * r(wannier_index,ir) * r(wannier_index,ir) )
										* exp( -alfa_new * r(wannier_index,ir) * 1.0/3.0 );
			}
		}
		
	}

	// 生成psir
	for(int i = 0; i < num_wannier; i++)
	{
		for(int ir = 0; ir < mesh_r; ir++)
		{
			psir(i,ir) = psi(i,ir) * r(i,ir);
		}
	}
	
	
	// 获得试探轨道
	for(int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
	{
		if(L[wannier_index] >= 0)
		{
			get_trial_orbitals_lm_k(wannier_index, L[wannier_index], m[wannier_index], ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
		}
		else
		{
			if(L[wannier_index] == -1 && m[wannier_index] == 0)
			{	
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> *tem_array = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs2 * tem_array[ig] + bs2 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array;
				
			}
			else if(L[wannier_index] == -1 && m[wannier_index] == 1)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> *tem_array = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs2 * tem_array[ig] - bs2 * trial_orbitals_k(wannier_index,ig);
				}	
				delete[] tem_array;
			}
			else if(L[wannier_index] == -2 && m[wannier_index] == 0)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] + bs2 * trial_orbitals_k(wannier_index,ig);
				}	
				delete[] tem_array_1;
				delete[] tem_array_2;
			}
			else if(L[wannier_index] == -2 && m[wannier_index] == 1)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] - bs2 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
			}			
			else if(L[wannier_index] == -2 && m[wannier_index] == 2)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs3 * tem_array[ig] + 2 * bs6 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array;
			}			
			else if(L[wannier_index] == -3 && m[wannier_index] == 0)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = 0.5*(tem_array_1[ig] + tem_array_2[ig] + tem_array_3[ig] + trial_orbitals_k(wannier_index,ig));
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
				
			}			
			else if(L[wannier_index] == -3 && m[wannier_index] == 1)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = 0.5*(tem_array_1[ig] + tem_array_2[ig] - tem_array_3[ig] - trial_orbitals_k(wannier_index,ig));
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -3 && m[wannier_index] == 2)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = 0.5*(tem_array_1[ig] - tem_array_2[ig] + tem_array_3[ig] - trial_orbitals_k(wannier_index,ig));
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -3 && m[wannier_index] == 3)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = 0.5*(tem_array_1[ig] - tem_array_2[ig] - tem_array_3[ig] + trial_orbitals_k(wannier_index,ig));
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -4 && m[wannier_index] == 0)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] + bs2 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
			}			
			else if(L[wannier_index] == -4 && m[wannier_index] == 1)
			{	
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] - bs2 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
			}			
			else if(L[wannier_index] == -4 && m[wannier_index] == 2)
			{	
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs3 * tem_array_1[ig] - 2 * bs6 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
			}			
			else if(L[wannier_index] == -4 && m[wannier_index] == 3)
			{
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs2 * tem_array_1[ig] + bs2 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
			}			
			else if(L[wannier_index] == -4 && m[wannier_index] == 4)
			{	
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = -1.0 * bs2 * tem_array_1[ig] + bs2 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
			}			
			else if(L[wannier_index] == -5 && m[wannier_index] == 0)
			{	
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig] - bs12 * tem_array_3[ig] + 0.5 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -5 && m[wannier_index] == 1)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig] - bs12 * tem_array_3[ig] + 0.5 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -5 && m[wannier_index] == 2)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig] - bs12 * tem_array_3[ig] - 0.5 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -5 && m[wannier_index] == 3)
			{	
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_3 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_3[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig] - bs12 * tem_array_3[ig] - 0.5 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
				delete[] tem_array_3;
			}			
			else if(L[wannier_index] == -5 && m[wannier_index] == 4)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig] + bs3 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
			}			
			else if(L[wannier_index] == -5 && m[wannier_index] == 5)
			{
				get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_1 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_1[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				complex<double> * tem_array_2 = new complex<double>[npwx];
				for(int ig = 0; ig < npwx; ig++)
				{
					tem_array_2[ig] = trial_orbitals_k(wannier_index,ig);
				}
				get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr,r,psir,mesh_r,gk,npw,trial_orbitals_k);
				for(int ig = 0; ig < npwx; ig++)
				{
					trial_orbitals_k(wannier_index,ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig] + bs3 * trial_orbitals_k(wannier_index,ig);
				}
				delete[] tem_array_1;
				delete[] tem_array_2;
			}	
		}
	}

	
	
}

// 注意这里轨道的L值必须是大于等于0的
void toWannier90::get_trial_orbitals_lm_k(const int wannier_index, const int orbital_L, const int orbital_m, matrix &ylm, 
										matrix &dr, matrix &r, matrix &psir, const int mesh_r, 
										Vector3<double> *gk, const int npw, ComplexMatrix &trial_orbitals_k)
{
	//计算径向函数在某个k点下倒空间的投影
	double *psik = new double[npw];
	double *psir_tem = new double[mesh_r];
	double *r_tem = new double[mesh_r];
	double *dr_tem = new double[mesh_r];
	double *psik_tem = new double[NQX];    //径向函数在固定k空间的投影（临时使用的数组）
	ZEROS(psir_tem,mesh_r);
	ZEROS(r_tem,mesh_r);
	ZEROS(dr_tem,mesh_r);
	
	for(int ir = 0; ir < mesh_r; ir++)
	{
		psir_tem[ir] = psir(wannier_index,ir);
		r_tem[ir] = r(wannier_index,ir);
		dr_tem[ir] = dr(wannier_index,ir);
	}
	
	toWannier90::integral(mesh_r,psir_tem,r_tem,dr_tem,orbital_L,psik_tem);
	
	// 从NQX个G点中插值法获得npw个G点的值
	for(int ig = 0; ig < npw; ig++)
	{
		psik[ig] = Mathzone::Polynomial_Interpolation(psik_tem, NQX, DQ, gk[ig].norm() * ucell.tpiba);
	}
	
	
	// 2.计算与原点选择（即轨道中心）而产生的相位在平面波基组下	
	complex<double> *sk = new complex<double>[npw];
	for(int ig = 0; ig < npw; ig++)
	{
		const double arg = ( gk[ig] * R_centre[wannier_index] ) * TWO_PI;
		sk[ig] = complex <double> ( cos(arg),  -sin(arg) );
	}
	
	// 3.生成 wannier_ylm
	double *wannier_ylm = new double[npw];
	for(int ig = 0; ig < npw; ig++)
	{
		int index = orbital_L * orbital_L + orbital_m;
		if(index == 2 || index == 3 || index == 5 || index == 6 || index == 14 || index == 15)
		{
			wannier_ylm[ig] = -1 * ylm(index,ig);
		}
		else
		{
			wannier_ylm[ig] = ylm(index,ig);
		}
	}
	
	// 4.计算最终试探轨道在某个k点下平面波基组的投影
	complex<double> lphase = pow(NEG_IMAG_UNIT, orbital_L);
	for(int ig = 0; ig < wf.npwx; ig++)
	{
		if(ig < npw)
		{
			trial_orbitals_k(wannier_index,ig) = lphase * sk[ig] * wannier_ylm[ig] * psik[ig];
		}
		else trial_orbitals_k(wannier_index,ig) = complex<double>(0.0,0.0);
	}
	
	
	// 5.归一化
	complex<double> anorm(0.0,0.0);
	for(int ig = 0; ig < wf.npwx; ig++)
	{
		anorm = anorm + conj(trial_orbitals_k(wannier_index,ig)) * trial_orbitals_k(wannier_index,ig);
	}
	
	complex<double> anorm_tem(0.0,0.0);
	MPI_Allreduce(&anorm , &anorm_tem , 1, MPI_DOUBLE_COMPLEX , MPI_SUM , POOL_WORLD);
	
	for(int ig = 0; ig < wf.npwx; ig++)
	{
		trial_orbitals_k(wannier_index,ig) = trial_orbitals_k(wannier_index,ig) / sqrt(anorm_tem);
	}
	
	delete[] psik;
	delete[] psir_tem;
	delete[] r_tem;
	delete[] dr_tem;
	delete[] psik_tem;
	delete[] sk;
	delete[] wannier_ylm;
	
	return;
	
}


void toWannier90::integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double* table)
{
	const double pref = FOUR_PI / sqrt(ucell.omega);
	
	double *inner_part = new double[meshr];
	for(int ir=0; ir<meshr; ir++)
	{
		inner_part[ir] = psir[ir] * psir[ir];
	}
	
	double unit = 0.0;
	Mathzone::Simpson_Integral(meshr, inner_part, rab, unit);
	delete[] inner_part;

	double *aux = new double[meshr];
	double *vchi = new double[meshr];
	for (int iq=0; iq<NQX; iq++)
	{
		const double q = DQ * iq;
		Mathzone::Spherical_Bessel(meshr, r, q, l, aux);
		for (int ir = 0;ir < meshr;ir++)
		{
			vchi[ir] = psir[ir] * aux[ir] * r[ir];
		}
		
		double vqint = 0.0;
		Mathzone::Simpson_Integral(meshr, vchi, rab, vqint);

		table[iq] =  vqint * pref;
	}
	delete[] aux;
	delete[] vchi;
	return;
}


void toWannier90::ToRealSpace(const int &ik, const int &ib, const ComplexMatrix *evc, complex<double> *psir, const Vector3<double> G)
{
	// (1) set value
	complex<double> *phase = UFFT.porter;
    ZEROS( psir, pw.nrxx );
	ZEROS( phase, pw.nrxx);


    for (int ig = 0; ig < kv.ngk[ik]; ig++)
    {
        psir[ pw.ig2fftw[ wf.igk(ik,ig) ] ] = evc[ik](ib, ig);
    }
	
	// get the phase value in realspace
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		if (pw.gdirect[ig] == G)
		{
			phase[ pw.ig2fftw[ig] ] = complex<double>(1.0,0.0);
			break;
		}
	}
	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psir, 1);
	pw.FFT_wfc.FFT3D(phase, 1);
	

	
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psir[ir] = psir[ir] * phase[ir];
	}
    return;
}

complex<double> toWannier90::unkdotb(const complex<double> *psir, const int ikb, const int bandindex, const ComplexMatrix *wfc_pw)
{
	complex<double> result(0.0,0.0);
	int knumber = kv.ngk[ikb];
	complex<double> *porter = UFFT.porter;
	ZEROS( porter, pw.nrxx);
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		porter[ir] = psir[ir];
	}
	pw.FFT_wfc.FFT3D( porter, -1);
	
	
	for (int ig = 0; ig < knumber; ig++)
	{
		result = result + conj( porter[ pw.ig2fftw[wf.igk(ikb, ig)] ] ) * wfc_pw[ikb](bandindex,ig);	
		
	}
	return result;
}

complex<double> toWannier90::unkdotkb(const int &ik, const int &ikb, const int &iband_L, const int &iband_R, const Vector3<double> G, const ComplexMatrix *wfc_pw)
{
	// (1) set value
	complex<double> result(0.0,0.0);
	complex<double> *psir = new complex<double>[pw.nrxx];
	complex<double> *phase = UFFT.porter;
    ZEROS( psir, pw.nrxx );
	ZEROS( phase, pw.nrxx);


    for (int ig = 0; ig < kv.ngk[ik]; ig++)
    {
        psir[ pw.ig2fftw[ wf.igk(ik,ig) ] ] = wfc_pw[ik](iband_L, ig);
    }
	
	// get the phase value in realspace
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		if (pw.gdirect[ig] == G)
		{
			phase[ pw.ig2fftw[ig] ] = complex<double>(1.0,0.0);
			break;
		}
	}
	
	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psir, 1);
	pw.FFT_wfc.FFT3D(phase, 1);
		
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psir[ir] = psir[ir] * phase[ir];
	}

	pw.FFT_wfc.FFT3D( psir, -1);
	
	complex<double> result_tem(0.0,0.0);
	
	for (int ig = 0; ig < kv.ngk[ikb]; ig++)
	{
		result_tem = result_tem + conj( psir[ pw.ig2fftw[wf.igk(ikb, ig)] ] ) * wfc_pw[ikb](iband_R,ig);	
		
	}
	
	MPI_Allreduce(&result_tem , &result , 1, MPI_DOUBLE_COMPLEX , MPI_SUM , POOL_WORLD);	
	
	delete[] psir;	
	return result;	
	
}

complex<double> toWannier90::gamma_only_cal(const int &ib_L, const int &ib_R, const ComplexMatrix *wfc_pw, const Vector3<double> G)
{
	complex<double> *phase = new complex<double>[pw.nrxx];
	complex<double> *psir = new complex<double>[pw.nrxx];
	complex<double> *psir_2 = new complex<double>[pw.nrxx];
	ZEROS( phase, pw.nrxx);
	ZEROS( psir, pw.nrxx);
	ZEROS( psir_2, pw.nrxx);

    for (int ig = 0; ig < kv.ngk[0]; ig++)
    {
        //psir[ pw.ig2fftw[ wf.igk(0,ig) ] ] = wfc_pw[0](ib_L, ig);
		psir[ pw.ig2fftw[ wf.igk(0,ig) ] ] = complex<double> ( abs(wfc_pw[0](ib_L, ig)), 0.0 );
    }
	
	// get the phase value in realspace
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		if (pw.gdirect[ig] == G)
		{
			phase[ pw.ig2fftw[ig] ] = complex<double>(1.0,0.0);
			break;
		}
	}
	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psir, 1);
	pw.FFT_wfc.FFT3D(phase, 1);
	
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psir_2[ir] = conj(psir[ir]) * phase[ir];
	}
	
		for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psir[ir] = psir[ir] * phase[ir];
	}
	
	pw.FFT_wfc.FFT3D( psir, -1);
	pw.FFT_wfc.FFT3D( psir_2, -1);
	
	complex<double> result(0.0,0.0);
	
	for (int ig = 0; ig < kv.ngk[0]; ig++)
	{
		//result = result + conj(psir_2[ pw.ig2fftw[wf.igk(0,ig)] ]) * wfc_pw[0](ib_R,ig) + psir[ pw.ig2fftw[ wf.igk(0,ig)] ] * conj(wfc_pw[0](ib_R,ig));
		//complex<double> tem = complex<double>( abs(wfc_pw[0](ib_R,ig)), 0.0 );
		result = result +  conj(psir[ pw.ig2fftw[ wf.igk(0,ig)] ]);// * tem;
	}
	
	delete[] phase;
	delete[] psir;
	delete[] psir_2;
	
	return result;
	
}

//使用lcao_in_pw方法将lcao基组转成pw基组
void toWannier90::lcao2pw_basis(const int ik, ComplexMatrix &orbital_in_G)
{
	this->table_local.create(ucell.ntype, ucell.nmax_total, NQX);
	Wavefunc_in_pw::make_table_q(ORB.orbital_file, this->table_local);
	Wavefunc_in_pw::produce_local_basis_in_pw(ik, orbital_in_G, this->table_local);
}

// 从lcao基组下产生pw基组的波函数周期部分unk的值，unk_inLcao[ik](ib,ig),ig的范围是kv.ngk[ik]
void toWannier90::getUnkFromLcao()
{
	complex<double>*** lcao_wfc_global = new complex<double>**[num_kpts];
	for(int ik = 0; ik < num_kpts; ik++)
	{
		lcao_wfc_global[ik] = new complex<double>*[NBANDS];
		for(int ib = 0; ib < NBANDS; ib++)
		{
			lcao_wfc_global[ik][ib] = new complex<double>[NLOCAL];
			ZEROS(lcao_wfc_global[ik][ib], NLOCAL);
		}
	}
	
	
	
	this->unk_inLcao = new ComplexMatrix[num_kpts];
	ComplexMatrix *orbital_in_G = new ComplexMatrix[num_kpts];

	for(int ik = 0; ik < num_kpts; ik++)
	{
		// 获取全局的lcao的波函数系数
		get_lcao_wfc_global_ik(lcao_wfc_global[ik],LOWF.WFC_K[ik]);
	
		int npw = kv.ngk[ik];
		unk_inLcao[ik].create(NBANDS,wf.npwx);
		orbital_in_G[ik].create(NLOCAL,npw);
		this->lcao2pw_basis(ik,orbital_in_G[ik]);
	
	}
	
	// 将lcao基组的unk转成pw基组下的unk
	for(int ik = 0; ik < num_kpts; ik++)
	{
		for(int ib = 0; ib < NBANDS; ib++)
		{
			for(int ig = 0; ig < kv.ngk[ik]; ig++)
			{
				for(int iw = 0; iw < NLOCAL; iw++)
				{
					unk_inLcao[ik](ib,ig) += orbital_in_G[ik](iw,ig)*lcao_wfc_global[ik][ib][iw];
				}
			}
		}
	}
	
	// 归一化
	for(int ik = 0; ik < num_kpts; ik++)
	{
		for(int ib = 0; ib < NBANDS; ib++)
		{
			complex<double> anorm(0.0,0.0);
			for(int ig = 0; ig < kv.ngk[ik]; ig++)
			{
				anorm = anorm + conj( unk_inLcao[ik](ib,ig) ) * unk_inLcao[ik](ib,ig);
			}
			
			complex<double> anorm_tem(0.0,0.0);
			MPI_Allreduce(&anorm , &anorm_tem , 1, MPI_DOUBLE_COMPLEX , MPI_SUM , POOL_WORLD);
			
			for(int ig = 0; ig < kv.ngk[ik]; ig++)
			{
				unk_inLcao[ik](ib,ig) = unk_inLcao[ik](ib,ig) / sqrt(anorm_tem);
			}
			
		}
	}
	
	
	for(int ik = 0; ik < kv.nkstot; ik++)
	{
		for(int ib = 0; ib < NBANDS; ib++)
		{
			delete[] lcao_wfc_global[ik][ib];
		}
		delete[] lcao_wfc_global[ik];
	}
	delete[] lcao_wfc_global;
	
	delete[] orbital_in_G;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	return;
}

// 获取全局的lcao的波函数系数
void toWannier90::get_lcao_wfc_global_ik(complex<double> **ctot, complex<double> **cc)
{
	complex<double>* ctot_send = new complex<double>[NBANDS*NLOCAL];

	MPI_Status status;

	for (int i=0; i<DSIZE; i++)
	{
		if (DRANK==0)
		{
			if (i==0)
			{
				// get the wave functions from 'ctot',
				// save them in the matrix 'c'.
				for (int iw=0; iw<NLOCAL; iw++)
				{
					const int mu_local = GridT.trace_lo[iw];
					if (mu_local >= 0)
					{
						for (int ib=0; ib<NBANDS; ib++)
						{
							//ctot[ib][iw] = cc[ib][mu_local];
							ctot_send[ib*NLOCAL+iw] = cc[ib][mu_local];
						}
					}
				}
			}
			else
			{
				int tag;
				// receive lgd2
				int lgd2 = 0;
				tag = i * 3;
				MPI_Recv(&lgd2, 1, MPI_INT, i, tag, DIAG_WORLD, &status);
				if(lgd2==0)
				{

				}
				else
				{
					// receive trace_lo2
					tag = i * 3 + 1;
					int* trace_lo2 = new int[NLOCAL];
					MPI_Recv(trace_lo2, NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

					// receive crecv
					complex<double>* crecv = new complex<double>[NBANDS*lgd2];
					ZEROS(crecv, NBANDS*lgd2);
					tag = i * 3 + 2;
					MPI_Recv(crecv,NBANDS*lgd2,mpicomplex,i,tag,DIAG_WORLD, &status);
				
					for (int ib=0; ib<NBANDS; ib++)
					{
						for (int iw=0; iw<NLOCAL; iw++)
						{
							const int mu_local = trace_lo2[iw];
							if (mu_local>=0)
							{
								//ctot[ib][iw] = crecv[mu_local*NBANDS+ib];
								ctot_send[ib*NLOCAL+iw] = crecv[mu_local*NBANDS+ib];
							}
						}
					}
				
					delete[] crecv;
					delete[] trace_lo2;
				}
			}
		}// end DRANK=0
		else if ( i == DRANK)
		{
			int tag;

			// send GridT.lgd
			tag = DRANK * 3;
			MPI_Send(&GridT.lgd, 1, MPI_INT, 0, tag, DIAG_WORLD);

			if(GridT.lgd != 0)
			{
				// send trace_lo
				tag = DRANK * 3 + 1;
				MPI_Send(GridT.trace_lo, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

				// send cc
				complex<double>* csend = new complex<double>[NBANDS*GridT.lgd];
				ZEROS(csend, NBANDS*GridT.lgd);

				for (int ib=0; ib<NBANDS; ib++)
				{
					for (int mu=0; mu<GridT.lgd; mu++)
					{
						csend[mu*NBANDS+ib] = cc[ib][mu];
					}
				}
			
				tag = DRANK * 3 + 2;
				MPI_Send(csend, NBANDS*GridT.lgd, mpicomplex, 0, tag, DIAG_WORLD);

			

				delete[] csend;

			}
		}// end i==DRANK
		MPI_Barrier(DIAG_WORLD);
	}

	MPI_Bcast(ctot_send,NBANDS*NLOCAL,mpicomplex,0,DIAG_WORLD);

	for(int ib = 0; ib < NBANDS; ib++)
	{
		for(int iw = 0; iw < NLOCAL; iw++)
		{
			ctot[ib][iw] = ctot_send[ib*NLOCAL+iw];
		}
	}

	delete[] ctot_send;

	return;
}




