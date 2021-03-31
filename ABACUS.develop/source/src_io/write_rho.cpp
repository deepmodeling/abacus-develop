#include "src_pw/charge.h"
#include "src_pw/global.h"

void Charge::write_rho(
	const double* rho_save, 
	const int &is, 
	const int &iter, 
	const string &fn, 
	const int &precision, 
	const bool for_plot)
{
    TITLE("Charge","write_rho");

    if (out_charge==0) 
	{
		return;
	}
	else if(iter % out_charge != 0) 
	{
		return; // mohan add 2010-05-22
	}
	
	time_t start, end;
	ofstream ofs;
	
	if(MY_RANK==0)
	{
		start = time(NULL);
    	
		ofs.open(fn.c_str());
    	if (!ofs)
    	{
        	WARNING("Charge::write_rho","Can't create Charge File!");
    	}	

		//ofs_running << "\n Output charge file." << endl;

		ofs << ucell.latName << endl;//1
		ofs << " " << ucell.lat0 * 0.529177 << endl;
		ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

		for(int it=0; it<ucell.ntype; it++)
		{
			ofs << " " << ucell.atoms[it].label;
		}
		ofs << endl;
		for(int it=0; it<ucell.ntype; it++)
		{
			ofs << " " << ucell.atoms[it].na;
		}
		ofs << endl;
		ofs << "Direct" << endl;

		for(int it=0; it<ucell.ntype; it++)
		{
			for(int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				ofs << " " << ucell.atoms[it].taud[ia].x
					<< " " << ucell.atoms[it].taud[ia].y
					<< " " << ucell.atoms[it].taud[ia].z << endl;
			}
		}


		if(for_plot)
		{

		}
		else
		{
			ofs << "\n  " << NSPIN;
			if(NSPIN==1 || NSPIN == 4)
			{
				ofs << "\n " << en.ef << " (fermi energy)";
			}
			else if(NSPIN==2)
			{
				if(is==0)ofs << "\n " << en.ef_up << " (fermi energy for spin=1)"; 
				else if(is==1)ofs << "\n " << en.ef_dw << " (fermi energy for spin=2)";
			}
			else
			{
				WARNING_QUIT("write_rho","check nspin!");
			}
		}
		ofs << "\n  " << pw.ncx << " " << pw.ncy << " " << pw.ncz << endl;

		ofs << setprecision(precision);
		ofs << scientific;

	}

	
#ifndef __MPI
	int count=0;
	for(int k=0; k<pw.ncz; k++)
	{
		for(int j=0; j<pw.ncy; j++)
		{
			for(int i=0; i<pw.ncx; i++)
			{
				if(count%8==0) ofs << "\n";
				ofs << " " << rho_save[i*pw.ncy*pw.ncz + j*pw.ncz + k];
				++count;
			}
		}
	}
#else
//	for(int ir=0; ir<pw.nrxx; ir++) chr.rho[0][ir]=1; // for testing
//	ofs_running << "\n RANK_IN_POOL = " << RANK_IN_POOL;
	
	// only do in the first pool.
	if(MY_POOL==0)
	{
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
//			ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
		}

		
		int count=0;
		int nxy = pw.ncx * pw.ncy;
		double* zpiece = new double[nxy];

		// save the rho one z by one z.
		for(int iz=0; iz<pw.ncz; iz++)
		{
			//	cout << "\n iz=" << iz << endl;
			// tag must be different for different iz.
			ZEROS(zpiece, nxy);
			int tag = iz;
			MPI_Status ierror;

			// case 1: the first part of rho in processor 0.
			if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
			{
				for(int ir=0; ir<nxy; ir++)
				{
					// mohan change to rho_save on 2012-02-10
					// because this can make our next restart calculation lead
					// to the same dr2 as the one saved.
					zpiece[ir] = rho_save[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
				}
			}
			// case 2: > first part rho: send the rho to 
			// processor 0.
			else if(which_ip[iz] == RANK_IN_POOL )
			{
				for(int ir=0; ir<nxy; ir++)
				{
					//						zpiece[ir] = rho[is][ir*num_z[RANK_IN_POOL]+iz];
					zpiece[ir] = rho_save[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
				}
				MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
			}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
			else if(RANK_IN_POOL==0)
			{
				MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
				//					ofs_running << "\n Receieve First number = " << zpiece[0];
			}

			// write data	
			if(MY_RANK==0)
			{
				//	ofs << "\niz=" << iz;
				// mohan update 2011-03-30
				for(int iy=0; iy<pw.ncy; iy++)
				{
					for(int ix=0; ix<pw.ncx; ix++)
					{
						if(count%8==0) ofs << "\n";
						ofs << " " << zpiece[ix*pw.ncy+iy];
						++count;
					}
				}
			}
		}// end iz
		delete[] zpiece;
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(MY_RANK==0) 
	{
		end = time(NULL);
		OUT_TIME("write_rho",start,end);
		ofs.close();
	}

    return;
}
