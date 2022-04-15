#include "../src_pw/charge.h"
#include "../src_pw/global.h"
#include "../module_base/element_name.h"

void Charge::write_rho_cube(
	const double* rho_save, 
	const int &is, 
	const std::string &fn, 
	const int &precision) 
{
    ModuleBase::TITLE("Charge","write_rho_cube");

	time_t start, end;
	std::ofstream ofs;
	
	if(GlobalV::MY_RANK==0)
	{
		start = time(NULL);
    	
		ofs.open(fn.c_str());
    	if (!ofs)
    	{
        	ModuleBase::WARNING("Charge::write_rho","Can't create Charge File!");
    	}	

		ofs << "Cubefile created from ABACUS SCF calculation" << std::endl;
		ofs << "Contains the selected quantity on a FFT grid" << std::endl;

		ofs << GlobalC::ucell.nat << " 0.0 0.0 0.0 " << std::endl;
		double fac=GlobalC::ucell.lat0;
		ofs << GlobalC::pw.ncx 
			<< " " << fac*GlobalC::ucell.latvec.e11/double(GlobalC::pw.ncx) 
			<< " " << fac*GlobalC::ucell.latvec.e12/double(GlobalC::pw.ncx) 
			<< " " << fac*GlobalC::ucell.latvec.e13/double(GlobalC::pw.ncx) << std::endl;
		ofs << GlobalC::pw.ncy 
			<< " " << fac*GlobalC::ucell.latvec.e21/double(GlobalC::pw.ncy) 
			<< " " << fac*GlobalC::ucell.latvec.e22/double(GlobalC::pw.ncy) 
			<< " " << fac*GlobalC::ucell.latvec.e23/double(GlobalC::pw.ncy) << std::endl;
		ofs << GlobalC::pw.ncz 
			<< " " << fac*GlobalC::ucell.latvec.e31/double(GlobalC::pw.ncz) 
			<< " " << fac*GlobalC::ucell.latvec.e32/double(GlobalC::pw.ncz) 
			<< " " << fac*GlobalC::ucell.latvec.e33/double(GlobalC::pw.ncz) << std::endl;

		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				//convert from label to atomic number
				int z = 0;
				for(int j=0; j!=ModuleBase::element_name.size(); j++)
					if (GlobalC::ucell.atoms[it].label == ModuleBase::element_name[j])
					{
						z=j+1;
						break;
					}
				ofs << " " << z << " " << z
					<< " " << fac*GlobalC::ucell.atoms[it].taud[ia].x
					<< " " << fac*GlobalC::ucell.atoms[it].taud[ia].y
					<< " " << fac*GlobalC::ucell.atoms[it].taud[ia].z << std::endl;
			}
		}

//		ofs << "\n  " << GlobalC::pw.ncx << " " << GlobalC::pw.ncy << " " << GlobalC::pw.ncz << std::endl;

		ofs << std::setprecision(precision);
		ofs << scientific;

	}

	
#ifndef __MPI
	int count=0;
	for(int i=0; i<GlobalC::pw.ncx; i++)
	{
		for(int j=0; j<GlobalC::pw.ncy; j++)
		{
			for(int k=0; k<GlobalC::pw.ncz; k++)
			{
				if(count%6==0) ofs << "\n";
				ofs << " " << rho_save[i*GlobalC::pw.ncy*GlobalC::pw.ncz + j*GlobalC::pw.ncz + k];
				++count;
			}
		}
	}
#else
//	for(int ir=0; ir<GlobalC::pw.nrxx; ir++) chr.rho[0][ir]=1; // for testing
//	ofs_running << "\n RANK_IN_POOL = " << RANK_IN_POOL;
	
	// only do in the first pool.
	if(GlobalV::MY_POOL==0)
	{
		int nxyz = GlobalC::pw.ncx * GlobalC::pw.ncy * GlobalC::pw.ncz;
		double* wfc = new double[nxyz];
		ModuleBase::GlobalFunc::ZEROS(wfc, nxyz);

		// num_z: how many planes on processor 'ip'
    	int *num_z = new int[GlobalV::NPROC_IN_POOL];
    	ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    	for (int iz=0;iz<GlobalC::pw.nbz;iz++)
    	{
        	int ip = iz % GlobalV::NPROC_IN_POOL;
        	num_z[ip] += GlobalC::pw.bz;
    	}	

		// start_z: start position of z in 
		// processor ip.
    	int *start_z = new int[GlobalV::NPROC_IN_POOL];
    	ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    	for (int ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    	{
        	start_z[ip] = start_z[ip-1]+num_z[ip-1];
    	}	
		delete[] num_z;

		// which_ip: found iz belongs to which ip.
		int *which_ip = new int[GlobalC::pw.ncz];
		ModuleBase::GlobalFunc::ZEROS(which_ip, GlobalC::pw.ncz);
		for(int iz=0; iz<GlobalC::pw.ncz; iz++)
		{
			for(int ip=0; ip<GlobalV::NPROC_IN_POOL; ip++)
			{
				if(iz>=start_z[GlobalV::NPROC_IN_POOL-1]) 
				{
					which_ip[iz] = GlobalV::NPROC_IN_POOL-1;
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

		int nxy = GlobalC::pw.ncx * GlobalC::pw.ncy;
		double* zpiece = new double[nxy];

		// save the rho one z by one z.
		for(int iz=0; iz<GlobalC::pw.ncz; iz++)
		{
			//	std::cout << "\n iz=" << iz << std::endl;
			// tag must be different for different iz.
			ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
			int tag = iz;
			MPI_Status ierror;

			// case 1: the first part of rho in processor 0.
			if(which_ip[iz] == 0 && GlobalV::RANK_IN_POOL ==0)
			{
				for(int ir=0; ir<nxy; ir++)
				{
					// mohan change to rho_save on 2012-02-10
					// because this can make our next restart calculation lead
					// to the same scf_thr as the one saved.
					zpiece[ir] = rho_save[ir*GlobalC::pw.nczp+iz-start_z[GlobalV::RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*GlobalC::pw.nczp+iz=" << ir*GlobalC::pw.nczp+iz;
				}
			}
			// case 2: > first part rho: send the rho to 
			// processor 0.
			else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
			{
				for(int ir=0; ir<nxy; ir++)
				{
					//						zpiece[ir] = rho[is][ir*num_z[RANK_IN_POOL]+iz];
					zpiece[ir] = rho_save[ir*GlobalC::pw.nczp+iz-start_z[GlobalV::RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*GlobalC::pw.nczp+iz=" << ir*GlobalC::pw.nczp+iz;
				}
				MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
			}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
			else if(GlobalV::RANK_IN_POOL==0)
			{
				MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
				//					ofs_running << "\n Receieve First number = " << zpiece[0];
			}

			if(GlobalV::MY_RANK==0)
			{
				for(int ir=0; ir<nxy; ir++)
				{
					wfc[ir+iz*nxy]=zpiece[ir];
				}
			}
		}// end iz
		delete[] zpiece;

		// write data	
		if(GlobalV::MY_RANK==0)
		{
			//	ofs << "\niz=" << iz;
			// mohan update 2011-03-30
			//int count=0;

			for(int ix=0; ix<GlobalC::pw.ncx; ix++)
			{
				for(int iy=0; iy<GlobalC::pw.ncy; iy++)
				{
					for (int iz=0; iz<GlobalC::pw.ncz; iz++)
					{
						ofs << " " << wfc[iz*GlobalC::pw.ncx*GlobalC::pw.ncy+ix*GlobalC::pw.ncy+iy];
						if(iz%6==5 && iz!=GlobalC::pw.ncz-1) ofs << "\n";
						//++count;
					}
					ofs << "\n";
				}
			}
		}
		delete[] wfc;
	}//end mypool=0
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(GlobalV::MY_RANK==0) 
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_rho",start,end);
		ofs.close();
	}

    return;
}
