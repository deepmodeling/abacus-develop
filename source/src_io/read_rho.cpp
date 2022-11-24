#include "../src_pw/charge.h"
#include "../src_pw/global.h"

bool Charge::read_rho(const int &is, const std::string &fn, double* rho) //add by dwan
{
    ModuleBase::TITLE("Charge","read_rho");
    std::ifstream ifs(fn.c_str());
    if (!ifs) 
	{
		ModuleBase::WARNING_QUIT("read_rho","!!! Couldn't find the charge file !!! The default directory \n of SPIN1_CHG is OUT.suffix, or you must set read_file_dir \n to a specific directory. ");
	}
	else
	{
    	GlobalV::ofs_running << " Find the file, try to read charge from file." << std::endl;
	}

	bool quit=false;

    std::string name;
	ifs >> name;
    
	// check lattice constant, unit is Angstrom
	ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.lat0 * 0.529177,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e11,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e12,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e13,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e21,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e22,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e23,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e31,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e32,quit);
    ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e33,quit);

	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		ModuleBase::CHECK_STRING(ifs,GlobalC::ucell.atoms[it].label,quit);
	}

	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].na,quit);
	}

	std::string coordinate;
	ifs >> coordinate;

	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].x,quit);
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].y,quit);
			ModuleBase::CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].z,quit);
		}
	}

	if(GlobalV::NSPIN != 4) ModuleBase::CHECK_INT(ifs, GlobalV::NSPIN);
	else
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalV::PRENSPIN);
	}
	if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::en.ef);
		GlobalV::ofs_running << " read in fermi energy = " << GlobalC::en.ef << std::endl;
	}
	else if(GlobalV::NSPIN == 2)
	{
		if(is==0)ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::en.ef_up);
		else if(is==1)ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::en.ef_dw);
	}
	else 
	{
		ModuleBase::WARNING_QUIT("read_rho","check nspin!");
	}
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->nx);	
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->ny);	
	ModuleBase::CHECK_INT(ifs, GlobalC::rhopw->nz);	

#ifndef __MPI
	GlobalV::ofs_running << " Read SPIN = " << is+1 << " charge now." << std::endl;
	for(int k=0; k<GlobalC::rhopw->nz; k++)
	{
		// consistent with the write_rho, something is
		// wrong.... but it works now.
		for(int j=0; j<GlobalC::rhopw->ny; j++)
		{
			for(int i=0; i<GlobalC::rhopw->nx; i++)
			{
				ifs >> rho[i*GlobalC::rhopw->ny*GlobalC::rhopw->nz + j*GlobalC::rhopw->nz +k];
			}
		}
	}
#else
	
	const int nxy = GlobalC::rhopw->nx * GlobalC::rhopw->ny;
	double *zpiece = new double[nxy];
	for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
	{
		ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
		if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0))
		{
			//				GlobalV::ofs_running << " Read charge density iz=" << iz << std::endl;
			for(int j=0; j<GlobalC::rhopw->ny; j++)
			{
				for(int i=0; i<GlobalC::rhopw->nx; i++)
				{
					ifs >> zpiece[ i*GlobalC::rhopw->ny + j ];
				}
			}
		}
		GlobalC::Pgrid.zpiece_to_all(zpiece, iz, rho);
	}// iz
	delete[] zpiece;
#endif

    if(GlobalV::MY_RANK==0||(GlobalV::ESOLVER_TYPE == "sdft"&&GlobalV::RANK_IN_STOGROUP==0)) ifs.close();
    return true;
}
