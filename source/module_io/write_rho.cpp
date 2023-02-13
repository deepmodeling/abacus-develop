#include "module_io/rho_io.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/element_name.h"

void ModuleIO::write_rho(
	const double* rho_save, 
	const int &is, 
	const int &iter, 
	const std::string &fn, 
	const int &precision, 
	const bool for_plot)
{
	ModuleBase::TITLE("ModuleIO","write_rho");

	if (GlobalV::out_chg==0) 
	{
		return;
	}
	else if(iter % GlobalV::out_chg != 0) 
	{
		return; // mohan add 2010-05-22
	}
	
	time_t start, end;
	std::ofstream ofs;
	std::ofstream ofs_cube;

	std::string fn_cube = fn;
	fn_cube.append(".cube");
	
	if(GlobalV::MY_RANK==0)
	{
		start = time(NULL);
    	
		ofs.open(fn.c_str());
		ofs_cube.open(fn_cube.c_str());
		
    		if (!ofs || !ofs_cube)
    		{
        			ModuleBase::WARNING("ModuleIO::write_rho","Can't create Charge File!");
    		}	

		//GlobalV::ofs_running << "\n Output charge file." << std::endl;

		ofs << GlobalC::ucell.latName << std::endl;//1
		ofs << " " << GlobalC::ucell.lat0 * 0.529177 << std::endl;
		ofs << " " << GlobalC::ucell.latvec.e11 << " " << GlobalC::ucell.latvec.e12 << " " << GlobalC::ucell.latvec.e13 << std::endl;
		ofs << " " << GlobalC::ucell.latvec.e21 << " " << GlobalC::ucell.latvec.e22 << " " << GlobalC::ucell.latvec.e23 << std::endl;
		ofs << " " << GlobalC::ucell.latvec.e31 << " " << GlobalC::ucell.latvec.e32 << " " << GlobalC::ucell.latvec.e33 << std::endl;

		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			ofs << " " << GlobalC::ucell.atoms[it].label;
		}
		ofs << std::endl;
		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			ofs << " " << GlobalC::ucell.atoms[it].na;
		}
		ofs << std::endl;
		ofs << "Direct" << std::endl;

		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				ofs << " " << GlobalC::ucell.atoms[it].taud[ia].x
					<< " " << GlobalC::ucell.atoms[it].taud[ia].y
					<< " " << GlobalC::ucell.atoms[it].taud[ia].z << std::endl;
			}
		}


		if(for_plot)
		{

		}
		else
		{
			ofs << "\n  " << GlobalV::NSPIN;
			if(GlobalV::NSPIN==1 || GlobalV::NSPIN == 4)
			{
				ofs << "\n " << GlobalC::en.ef << " (fermi energy)";
			}
			else if(GlobalV::NSPIN==2)
			{
				if (GlobalV::TWO_EFERMI)
				{
					if(is==0)ofs << "\n " << GlobalC::en.ef_up << " (fermi energy for spin=1)"; 
					else if(is==1)ofs << "\n " << GlobalC::en.ef_dw << " (fermi energy for spin=2)";
				}
				else
				{
					ofs << "\n " << GlobalC::en.ef << " (fermi energy)";
				}
			}
			else
			{
				ModuleBase::WARNING_QUIT("write_rho","check nspin!");
			}
		}
		ofs << "\n  " << GlobalC::rhopw->nx << " " << GlobalC::rhopw->ny << " " << GlobalC::rhopw->nz << std::endl;
		ofs << std::setprecision(precision);
		ofs << scientific;

		/// output header for cube file
		ofs_cube << "Cubefile created from ABACUS SCF calculation" << std::endl;
		ofs_cube << "Contains the selected quantity on a FFT grid" << std::endl;
		ofs_cube << GlobalC::ucell.nat << " 0.0 0.0 0.0 " << std::endl;
		double fac=GlobalC::ucell.lat0;
		ofs_cube << GlobalC::rhopw->nx 
			<< " " << fac*GlobalC::ucell.latvec.e11/double(GlobalC::rhopw->nx) 
			<< " " << fac*GlobalC::ucell.latvec.e12/double(GlobalC::rhopw->nx) 
			<< " " << fac*GlobalC::ucell.latvec.e13/double(GlobalC::rhopw->nx) << std::endl;
		ofs_cube << GlobalC::rhopw->ny 
			<< " " << fac*GlobalC::ucell.latvec.e21/double(GlobalC::rhopw->ny) 
			<< " " << fac*GlobalC::ucell.latvec.e22/double(GlobalC::rhopw->ny) 
			<< " " << fac*GlobalC::ucell.latvec.e23/double(GlobalC::rhopw->ny) << std::endl;
		ofs_cube << GlobalC::rhopw->nz 
			<< " " << fac*GlobalC::ucell.latvec.e31/double(GlobalC::rhopw->nz) 
			<< " " << fac*GlobalC::ucell.latvec.e32/double(GlobalC::rhopw->nz) 
			<< " " << fac*GlobalC::ucell.latvec.e33/double(GlobalC::rhopw->nz) << std::endl;

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
				ofs_cube << " " << z << " " << z
					<< " " << fac*GlobalC::ucell.atoms[it].tau[ia].x
					<< " " << fac*GlobalC::ucell.atoms[it].tau[ia].y
					<< " " << fac*GlobalC::ucell.atoms[it].tau[ia].z << std::endl;
			}
		}
		ofs_cube << std::setprecision(precision);
		ofs_cube << scientific;
		/// cube file
	}
	
#ifndef __MPI
	int count=0;
	for(int k=0; k<GlobalC::rhopw->nz; k++)
	{
		for(int j=0; j<GlobalC::rhopw->ny; j++)
		{
			for(int i=0; i<GlobalC::rhopw->nx; i++)
			{
				if(count%8==0) ofs << "\n";
				ofs << " " << rho_save[i*GlobalC::rhopw->ny*GlobalC::rhopw->nz + j*GlobalC::rhopw->nz + k];
				++count;
			}
		}
	}
	/// for cube file
	int count_cube=0;
	for(int i=0; i<GlobalC::rhopw->nx; i++)
	{
		for(int j=0; j<GlobalC::rhopw->ny; j++)
		{
			for(int k=0; k<GlobalC::rhopw->nz; k++)
			{
				if(count_cube%6==0) ofs_cube << "\n";
				ofs_cube << " " << rho_save[i*GlobalC::rhopw->ny*GlobalC::rhopw->nz + j*GlobalC::rhopw->nz + k];
				++count_cube;
			}
		}
	}
	/// for cube file
#else
//	for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++) chr.rho[0][ir]=1; // for testing
//	GlobalV::ofs_running << "\n GlobalV::RANK_IN_POOL = " << GlobalV::RANK_IN_POOL;
	
	// only do in the first pool.
	if(GlobalV::MY_POOL==0)
	{
		/// for cube file
		int nxyz = GlobalC::rhopw->nx * GlobalC::rhopw->ny * GlobalC::rhopw->nz;
		double* chg_cube = new double[nxyz];
		ModuleBase::GlobalFunc::ZEROS(chg_cube, nxyz);
		/// for cube file
	
		// num_z: how many planes on processor 'ip'
    		int *num_z = new int[GlobalV::NPROC_IN_POOL];
    		ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    		for (int iz=0;iz<GlobalC::bigpw->nbz;iz++)
    		{
        		int ip = iz % GlobalV::NPROC_IN_POOL;
        		num_z[ip] += GlobalC::bigpw->bz;
    		}

		// start_z: start position of z in 
		// processor ip.
    		int *start_z = new int[GlobalV::NPROC_IN_POOL];
    		ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    		for (int ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    		{
        		start_z[ip] = start_z[ip-1]+num_z[ip-1];
    		}

		// which_ip: found iz belongs to which ip.
		int *which_ip = new int[GlobalC::rhopw->nz];
		ModuleBase::GlobalFunc::ZEROS(which_ip, GlobalC::rhopw->nz);
		for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
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
			// GlobalV::ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
		}

		
		int count=0;
		int nxy = GlobalC::rhopw->nx * GlobalC::rhopw->ny;
		double* zpiece = new double[nxy];

		// save the rho one z by one z.
		for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
		{
			// std::cout << "\n iz=" << iz << std::endl;
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
					zpiece[ir] = rho_save[ir*GlobalC::rhopw->nplane+iz-GlobalC::rhopw->startz_current];
					// GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*GlobalC::rhopw->nplane+iz=" << ir*GlobalC::rhopw->nplane+iz;
				}
			}
			// case 2: > first part rho: send the rho to 
			// processor 0.
			else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
			{
				for(int ir=0; ir<nxy; ir++)
				{
					// zpiece[ir] = rho[is][ir*num_z[GlobalV::RANK_IN_POOL]+iz];
					zpiece[ir] = rho_save[ir*GlobalC::rhopw->nplane+iz-GlobalC::rhopw->startz_current];
					// GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*GlobalC::rhopw->nplane+iz=" << ir*GlobalC::rhopw->nplane+iz;
				}
				MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
			}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
			else if(GlobalV::RANK_IN_POOL==0)
			{
				MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
				// GlobalV::ofs_running << "\n Receieve First number = " << zpiece[0];
			}

			// write data	
			if(GlobalV::MY_RANK==0)
			{
				// ofs << "\niz=" << iz;
				// mohan update 2011-03-30
				for(int iy=0; iy<GlobalC::rhopw->ny; iy++)
				{
					for(int ix=0; ix<GlobalC::rhopw->nx; ix++)
					{
						if(count%8==0) ofs << "\n";
						ofs << " " << zpiece[ix*GlobalC::rhopw->ny+iy];
						++count;
					}
				}
				/// for cube file
				for(int ir=0; ir<nxy; ir++)
				{
					chg_cube[ir+iz*nxy]=zpiece[ir];
				}
				/// for cube file
			}
		}// end iz
		delete[] zpiece;
		delete[] which_ip;
		delete[] num_z;
		delete[] start_z;
		// for cube file
		if(GlobalV::MY_RANK==0)
		{
			for(int ix=0; ix<GlobalC::rhopw->nx; ix++)
			{
				for(int iy=0; iy<GlobalC::rhopw->ny; iy++)
				{
					for (int iz=0; iz<GlobalC::rhopw->nz; iz++)
					{
						ofs_cube << " " << chg_cube[iz*GlobalC::rhopw->nx*GlobalC::rhopw->ny+ix*GlobalC::rhopw->ny+iy];
						if(iz%6==5 && iz!=GlobalC::rhopw->nz-1) ofs_cube << "\n";
					}
					ofs_cube << "\n";
				}
			}
		}
		delete[] chg_cube;
		/// for cube file
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(GlobalV::MY_RANK==0) 
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_rho",start,end);
		ofs.close();
		/// for cube file
		ofs_cube.close();
	}

    return;
}
