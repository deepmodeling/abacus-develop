#include "dmk_io.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"

// output DM_k
void ModuleIO::write_dmk(
    const K_Vectors& kv,
    const int& ik,
	const std::string &fn, 
	const int &precision,
	std::vector<ModuleBase::ComplexMatrix> &dm_k)
{
    ModuleBase::TITLE("ModuleIO","write_dm");

	ModuleBase::timer::tick("ModuleIO","write_dm");

	time_t start, end;
	std::ofstream ofs;

	if(GlobalV::MY_RANK==0)
	{
		start = time(NULL);

		ofs.open(fn.c_str());
		if (!ofs)
		{
			ModuleBase::WARNING("ModuleIO::write_dm","Can't create DENSITY MATRIX File!");
		}

        ofs << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y << " " << kv.kvec_d[ik].z << std::endl;
		//ofs << "\n " << ef << " (fermi energy)";

		ofs << "\n  " << GlobalV::NLOCAL << " " << GlobalV::NLOCAL << std::endl;

		ofs << std::setprecision(precision);
		ofs << std::scientific;

	}

    for(int i=0; i<GlobalV::NLOCAL; ++i)
    {
        for(int j=0; j<GlobalV::NLOCAL; ++j)
        {
            if(j%8==0) ofs << "\n";
            ofs << " " << dm_k[ik](i,j).real();
            //ofs << " " << DM[is][i][j];
        }
    }

	if(GlobalV::MY_RANK==0)
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_dm",start,end);
		ofs.close();
	}
	ModuleBase::timer::tick("ModuleIO","write_dm");

    return;
}


void ModuleIO::read_dmk(
    const K_Vectors& kv,
    const int& ik,
	const std::string &fn,
	std::vector<ModuleBase::ComplexMatrix> &dm_k)
{
    //weiqing modify 2023-8-7
    bool quit_abacus = false;

    std::ifstream ifs;
    if(GlobalV::MY_RANK==0)
    {
        ifs.open(fn.c_str());
        if (!ifs)
        {
            //xiaohui modify 2015-03-25
            //quit_mesia = true;
            quit_abacus = true;
        }
        else
        {
            // if the number is not match,
            // quit the program or not.
            bool quit=false;

            ModuleBase::CHECK_DOUBLE(ifs,kv.kvec_d[ik].x,quit);
            ModuleBase::CHECK_DOUBLE(ifs,kv.kvec_d[ik].y,quit);
            ModuleBase::CHECK_DOUBLE(ifs,kv.kvec_d[ik].z,quit);

            ModuleBase::CHECK_INT(ifs, GlobalV::NLOCAL);
            ModuleBase::CHECK_INT(ifs, GlobalV::NLOCAL);
        }// If file exist, read in data.
    } // Finish reading the first part of density matrix.

    for(int i=0; i<GlobalV::NLOCAL; ++i)
    {
        for(int j=0; j<GlobalV::NLOCAL; ++j)
        {
            ifs >> dm_k[ik](i,j);
        }
    }

    ifs.close();

    return;
}
