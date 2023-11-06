#include "module_base/element_name.h"
#include "module_io/cube_io.h"

void ModuleIO::write_cube(
#ifdef __MPI
	const int& bz,
	const int& nbz,
	const int& nplane,
	const int& startz_current,
#endif
	const double* data,
	const int& is,
	const int& nspin,
	const int& iter,
	const std::string& fn,
	const int& nx,
	const int& ny,
	const int& nz,
	const double& ef,
	const UnitCell* ucell,
	const int &precision,
	const int &out_fermi)
{
	ModuleBase::TITLE("ModuleIO","write_cube");

	time_t start, end;
	std::ofstream ofs_cube;
  
	if(GlobalV::MY_RANK==0)
	{
		start = time(NULL);

		ofs_cube.open(fn.c_str());
		
		if (!ofs_cube)
		{
				ModuleBase::WARNING("ModuleIO::write_cube","Can't create Output File!");
		}	

		/// output header for cube file
		ofs_cube << "Cubefile created from ABACUS SCF calculation. The inner loop is z index, followed by y index, x index in turn." << std::endl;
		// ofs_cube << "Contains the selected quantity on a FFT grid" << std::endl;
		ofs_cube << nspin << " (nspin) ";
		
		ofs_cube << std::fixed;
		ofs_cube << std::setprecision(6);
		if (out_fermi == 1)
		{
			if (GlobalV::TWO_EFERMI)
			{
				if(is==0)	ofs_cube << ef << " (fermi energy for spin=1, in Ry)" << std::endl;
				else if(is==1)	ofs_cube << ef << " (fermi energy for spin=2, in Ry)" << std::endl;
			}
			else
			{
				ofs_cube << ef << " (fermi energy, in Ry)" << std::endl;
			}
		}
		else
		{
			ofs_cube << std::endl;
		}
		
		ofs_cube << ucell->nat << " 0.0 0.0 0.0 " << std::endl;
		double fac=ucell->lat0;
		ofs_cube << nx 
			<< " " << fac*ucell->latvec.e11/double(nx) 
			<< " " << fac*ucell->latvec.e12/double(nx) 
			<< " " << fac*ucell->latvec.e13/double(nx) << std::endl;
		ofs_cube << ny 
			<< " " << fac*ucell->latvec.e21/double(ny) 
			<< " " << fac*ucell->latvec.e22/double(ny) 
			<< " " << fac*ucell->latvec.e23/double(ny) << std::endl;
		ofs_cube << nz 
			<< " " << fac*ucell->latvec.e31/double(nz) 
			<< " " << fac*ucell->latvec.e32/double(nz) 
			<< " " << fac*ucell->latvec.e33/double(nz) << std::endl;

		std::string element = "";
		for(int it=0; it<ucell->ntype; it++)
		{
			// erase the number in label, such as Fe1.
			element = ucell->atoms[it].label;
			std::string::iterator temp = element.begin();
			while (temp != element.end())
			{
				if ((*temp >= '1') && (*temp <= '9'))
				{
					temp = element.erase(temp);
				}
				else
				{
					temp++;
				}
			}

			for(int ia=0; ia<ucell->atoms[it].na; ia++)
			{
				//convert from label to atomic number
				int z = 0;
				for(int j=0; j!=ModuleBase::element_name.size(); j++)
				{
					if (element == ModuleBase::element_name[j])
					{
						z=j+1;
						break;
					}
				}
				ofs_cube << " " << z << " " << ucell->atoms[it].ncpp.zv
						 << " " << fac*ucell->atoms[it].tau[ia].x
						 << " " << fac*ucell->atoms[it].tau[ia].y
						 << " " << fac*ucell->atoms[it].tau[ia].z << std::endl;
			}
		}
		ofs_cube.unsetf(std::ostream::fixed);
		ofs_cube << std::setprecision(precision);
		ofs_cube << std::scientific;
	}

#ifdef __MPI
//	for(int ir=0; ir<rhopw->nrxx; ir++) chr.rho[0][ir]=1; // for testing
//	GlobalV::ofs_running << "\n GlobalV::RANK_IN_POOL = " << GlobalV::RANK_IN_POOL;
	
	// only do in the first pool.
	if(GlobalV::MY_POOL==0)
	{
		/// for cube file
		int nxyz = nx * ny * nz;
		double* data_cube = new double[nxyz];
		ModuleBase::GlobalFunc::ZEROS(data_cube, nxyz);
		/// for cube file
	
		// num_z: how many planes on processor 'ip'
    		int *num_z = new int[GlobalV::NPROC_IN_POOL];
    		ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    		for (int iz=0;iz<nbz;iz++)
    		{
        		int ip = iz % GlobalV::NPROC_IN_POOL;
        		num_z[ip] += bz;
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
		int *which_ip = new int[nz];
		ModuleBase::GlobalFunc::ZEROS(which_ip, nz);
		for(int iz=0; iz<nz; iz++)
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
		int nxy = nx * ny;
		double* zpiece = new double[nxy];

		// save the rho one z by one z.
		for(int iz=0; iz<nz; iz++)
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
					zpiece[ir] = data[ir*nplane+iz-startz_current];
					// GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*rhopw->nplane+iz=" << ir*rhopw->nplane+iz;
				}
			}
			// case 2: > first part rho: send the rho to 
			// processor 0.
			else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
			{
				for(int ir=0; ir<nxy; ir++)
				{
					// zpiece[ir] = rho[is][ir*num_z[GlobalV::RANK_IN_POOL]+iz];
					zpiece[ir] = data[ir*nplane+iz-startz_current];
					// GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*rhopw->nplane+iz=" << ir*rhopw->nplane+iz;
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

			if(GlobalV::MY_RANK==0)
			{
				/// for cube file
				for(int ir=0; ir<nxy; ir++)
				{
					data_cube[ir+iz*nxy]=zpiece[ir];
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
			for(int ix=0; ix<nx; ix++)
			{
				for(int iy=0; iy<ny; iy++)
				{
					for (int iz=0; iz<nz; iz++)
					{
						ofs_cube << " " << data_cube[iz*nx*ny+ix*ny+iy];
						if(iz%6==5 && iz!=nz-1) ofs_cube << "\n";
					}
					ofs_cube << "\n";
				}
			}
		}
		delete[] data_cube;
		/// for cube file
	}
	MPI_Barrier(MPI_COMM_WORLD);
#else
	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int k=0; k<nz; k++)
			{
				ofs_cube << " " << data[k*nx*ny+i*ny+j];
				// ++count_cube;
				if(k%6==5 && k!=nz-1) ofs_cube << "\n";
			}
			ofs_cube << "\n";
		}
	}
#endif

	if(GlobalV::MY_RANK==0) 
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_cube",start,end);
		// ofs.close();
		/// for cube file
		ofs_cube.close();
	}

    return;
}
