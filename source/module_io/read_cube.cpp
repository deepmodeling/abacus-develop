#include "module_io/cube_io.h"
#include "module_base/global_variable.h"

bool ModuleIO::read_cube(
#ifdef __MPI
    Parallel_Grid* Pgrid,
#endif
    const int& is,
    const int& nspin,
    const std::string& fn,
    double* data,
    const int& nx,
    const int& ny,
    const int& nz,
    double& ef,
    const UnitCell* ucell,
    int& prenspin)
{
    ModuleBase::TITLE("ModuleIO","read_cube");
    std::ifstream ifs(fn.c_str());
    if (!ifs) 
	{
		std::string tmp_warning_info = "!!! Couldn't find the charge file of ";
		tmp_warning_info += fn;
		GlobalV::ofs_running << tmp_warning_info << std::endl;
		return false;
	}
	else
	{
    	GlobalV::ofs_running << " Find the file, try to read charge from file." << std::endl;
	}

	bool quit=false;

	ifs.ignore(300, '\n'); // skip the header

	if(nspin != 4)
	{
		ModuleBase::CHECK_INT(ifs, nspin);
	}
	else
	{
		ifs >> prenspin;
	}
	ifs.ignore(150, ')');

	ifs >> ef;
	GlobalV::ofs_running << " read in fermi energy = " << ef << std::endl;

	ifs.ignore(150, '\n');

	ModuleBase::CHECK_INT(ifs,ucell->nat,quit);
	ifs.ignore(150, '\n');

	double fac=ucell->lat0;
	ModuleBase::CHECK_INT(ifs, nx);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e11/double(nx), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e12/double(nx), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e13/double(nx), quit);
	ModuleBase::CHECK_INT(ifs, ny);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e21/double(ny), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e22/double(ny), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e23/double(ny), quit);
	ModuleBase::CHECK_INT(ifs, nz);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e31/double(nz), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e32/double(nz), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*ucell->latvec.e33/double(nz), quit);

	int temp = 0;
	for(int it=0; it<ucell->ntype; it++)
	{
		for(int ia=0; ia<ucell->atoms[it].na; ia++)
		{
			ifs >> temp; // skip atomic number
			ifs >> temp; // skip Z valance
			// ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell->atoms[it].ncpp.zv,quit); // check Z valance
			ModuleBase::CHECK_DOUBLE(ifs,fac*ucell->atoms[it].tau[ia].x,quit);
			ModuleBase::CHECK_DOUBLE(ifs,fac*ucell->atoms[it].tau[ia].y,quit);
			ModuleBase::CHECK_DOUBLE(ifs,fac*ucell->atoms[it].tau[ia].z,quit);
		}
	}

#ifdef __MPI
	const int nxy = nx * ny;
	double *zpiece = nullptr;
	double **tempRho = nullptr;

	if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
	{
		tempRho = new double*[nz];
		for(int iz=0; iz<nz; iz++)
		{
			tempRho[iz] = new double[nxy];
			// ModuleBase::GlobalFunc::ZEROS(tempRho[iz], nxy);
		}
		for(int ix=0; ix<nx; ix++)
		{
			for(int iy=0; iy<ny; iy++)
			{
				for(int iz=0; iz<nz; iz++)
				{
					ifs >> tempRho[iz][ix*ny + iy];
				}
			}
		}
	}
	else
	{
		zpiece = new double[nxy];
		ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
	}

	for(int iz=0; iz<nz; iz++)
	{
		if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
		{
			zpiece = tempRho[iz];
		}
		Pgrid->zpiece_to_all(zpiece, iz, data);
	}// iz

	if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
	{
		for(int iz=0; iz<nz; iz++)
		{
			delete[] tempRho[iz];
		}
		delete[] tempRho;
	}
	else
	{
		delete[] zpiece;
	}
#else
	GlobalV::ofs_running << " Read SPIN = " << is+1 << " charge now." << std::endl;
	// consistent with the write_cube,
	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int k=0; k<nz; k++)
			{
				ifs >> data[k*nx*ny+i*ny+j];
			}
		}
	}
#endif

    if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0)) ifs.close();
    return true;
}
