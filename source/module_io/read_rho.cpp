#include "module_io/rho_io.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

bool ModuleIO::read_rho(const int &is, const std::string &fn, double* rho, int &prenspin) //add by dwan
{
    ModuleBase::TITLE("ModuleIO","read_rho");
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

	if(GlobalV::NSPIN != 4) ModuleBase::CHECK_INT(ifs, GlobalV::NSPIN);
	else
	{
		ifs >> prenspin;
	}
	ifs.ignore(150, ')');

	if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
	{
		ifs >> GlobalC::en.ef;
		GlobalV::ofs_running << " read in fermi energy = " << GlobalC::en.ef << std::endl;
	}
	else if(GlobalV::NSPIN == 2)
	{
		if(is==0)		ifs >> GlobalC::en.ef_up;
		else if(is==1)	ifs >> GlobalC::en.ef_dw;
	}
	else 
	{
		ModuleBase::WARNING_QUIT("read_rho","check nspin!");
	}
	ifs.ignore(150, '\n');

	ModuleBase::CHECK_INT(ifs,GlobalC::ucell.nat,quit);
	ifs.ignore(150, '\n');

	double fac=GlobalC::ucell.lat0;
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->nx);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e11/double(GlobalC::rhopw->nx), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e12/double(GlobalC::rhopw->nx), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e13/double(GlobalC::rhopw->nx), quit);
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->ny);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e21/double(GlobalC::rhopw->ny), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e22/double(GlobalC::rhopw->ny), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e23/double(GlobalC::rhopw->ny), quit);
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->nz);	
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e31/double(GlobalC::rhopw->nz), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e32/double(GlobalC::rhopw->nz), quit);
	ModuleBase::CHECK_DOUBLE(ifs, fac*GlobalC::ucell.latvec.e33/double(GlobalC::rhopw->nz), quit);

	int temp = 0;
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			ifs >> temp; // skip atomic number
			ifs >> temp; // skip Z valance
			// ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].ncpp.zv,quit); // check Z valance
			ModuleBase::CHECK_DOUBLE(ifs,fac*GlobalC::ucell.atoms[it].taud[ia].x,quit);
			ModuleBase::CHECK_DOUBLE(ifs,fac*GlobalC::ucell.atoms[it].taud[ia].y,quit);
			ModuleBase::CHECK_DOUBLE(ifs,fac*GlobalC::ucell.atoms[it].taud[ia].z,quit);
		}
	}

#ifndef __MPI
	GlobalV::ofs_running << " Read SPIN = " << is+1 << " charge now." << std::endl;
	// consistent with the write_rho,
	for(int i=0; i<GlobalC::rhopw->nx; i++)
	{
		for(int j=0; j<GlobalC::rhopw->ny; j++)
		{
			for(int k=0; k<GlobalC::rhopw->nz; k++)
			{
				ifs >> rho[k*GlobalC::rhopw->nx*GlobalC::rhopw->ny+i*GlobalC::rhopw->ny+j];
			}
		}
	}
#else
	
	const int nxy = GlobalC::rhopw->nx * GlobalC::rhopw->ny;
	double *zpiece = nullptr;
	double **tempRho = nullptr;

	if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
	{
		tempRho = new double*[GlobalC::rhopw->nz];
		for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
		{
			tempRho[iz] = new double[nxy];
			// ModuleBase::GlobalFunc::ZEROS(tempRho[iz], nxy);
		}
		for(int ix=0; ix<GlobalC::rhopw->nx; ix++)
		{
			for(int iy=0; iy<GlobalC::rhopw->ny; iy++)
			{
				for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
				{
					ifs >> tempRho[iz][ix*GlobalC::rhopw->ny + iy];
				}
			}
		}
	}
	else
	{
		zpiece = new double[nxy];
		ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
	}

	for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
	{
		if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
		{
			zpiece = tempRho[iz];
		}
		GlobalC::Pgrid.zpiece_to_all(zpiece, iz, rho);
	}// iz

	if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
	{
		for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
		{
			delete[] tempRho[iz];
		}
		delete[] tempRho;
	}
	else
	{
		delete[] zpiece;
	}
#endif

    if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0)) ifs.close();
    return true;
}
