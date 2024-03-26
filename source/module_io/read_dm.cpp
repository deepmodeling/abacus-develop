#include "module_io/dm_io.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"


void ModuleIO::read_dm(
#ifdef __MPI
	const int nnrg,
	const int* trace_lo,
#endif
    const bool gamma_only_local,
    const int nlocal,
    const int nspin,
	const int &is,
	const std::string &fn,
	double*** DM,
	double** DM_R,
	double& ef,
	const UnitCell* ucell)
{
    ModuleBase::TITLE("ModuleIO","read_dm");
    ModuleBase::timer::tick("ModuleIO","read_dm");

    GlobalV::ofs_running << "\n processor 0 is reading density matrix from file < " << fn << " > " << std::endl;
    //xiaohui modify 2015-03-25
    //bool quit_mesia = false;
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

            std::string name;
            ifs >> name;

            // check lattice constant, unit is Angstrom
            ModuleBase::CHECK_DOUBLE(ifs,ucell->lat0 * ModuleBase::BOHR_TO_A,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e11,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e12,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e13,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e21,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e22,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e23,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e31,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e32,quit);
            ModuleBase::CHECK_DOUBLE(ifs,ucell->latvec.e33,quit);

            for(int it=0; it<ucell->ntype; it++)
            {
                ModuleBase::CHECK_STRING(ifs,ucell->atoms[it].label,quit);
            }

            for(int it=0; it<ucell->ntype; it++)
            {
                ModuleBase::CHECK_DOUBLE(ifs,ucell->atoms[it].na,quit);
            }

            std::string coordinate;
            ifs >> coordinate;

            for(int it=0; it<ucell->ntype; it++)
            {
                for(int ia=0; ia<ucell->atoms[it].na; ia++)
                {
                    ModuleBase::CHECK_DOUBLE(ifs,ucell->atoms[it].taud[ia].x,quit);
                    ModuleBase::CHECK_DOUBLE(ifs,ucell->atoms[it].taud[ia].y,quit);
                    ModuleBase::CHECK_DOUBLE(ifs,ucell->atoms[it].taud[ia].z,quit);
                }
            }

            ModuleBase::CHECK_INT(ifs, nspin);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ef);
            ModuleBase::CHECK_INT(ifs, nlocal);
            ModuleBase::CHECK_INT(ifs, nlocal);
        }// If file exist, read in data.
    } // Finish reading the first part of density matrix.


#ifndef __MPI
    GlobalV::ofs_running << " Read SPIN = " << is+1 << " density matrix now." << std::endl;

    if(gamma_only_local)
    {
        for(int i=0; i<nlocal; ++i)
        {
            for(int j=0; j<nlocal; ++j)
            {
                ifs >> DM[is][i][j];
            }
        }
    }
    else
    {
    #ifdef __MPI
        ModuleBase::WARNING_QUIT("ModuleIO::read_dm","The nnrg should not be update");
        ModuleBase::CHECK_INT(ifs,nnrg);

        for(int i=0; i<nnrg; ++i)
        {
            ifs >> DM_R[is][i];
        }
    #endif
    }
#else

    // distribution of necessary data
    //xiaohui modify 2015-03-25
    //Parallel_Common::bcast_bool(quit_mesia);
    Parallel_Common::bcast_bool(quit_abacus);
    //xiaohui modify 2015-03-25
    //if(quit_mesia)
    if(quit_abacus)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::read_dm","Can not find the density matrix file.");
    }

    Parallel_Common::bcast_double(ef);

    if(gamma_only_local)
    {

        double *tmp = new double[nlocal];
        for(int i=0; i<nlocal; ++i)
        {
            //GlobalV::ofs_running << " i=" << i << std::endl;
            ModuleBase::GlobalFunc::ZEROS(tmp, nlocal);
            if(GlobalV::MY_RANK==0)
            {
                for(int j=0; j<nlocal; ++j)
                {
                    ifs >> tmp[j];
                }
            }
            Parallel_Common::bcast_double(tmp, nlocal);

            const int mu = trace_lo[i];
            if(mu >= 0)
            {   
                for(int j=0; j<nlocal; ++j)
                {
                    const int nu = trace_lo[j];
                    if(nu >= 0)
                    {
                        DM[is][mu][nu] = tmp[j];
                    }
                }
            }
        }// i
        delete[] tmp;
    }
    else
    {
        ModuleBase::WARNING_QUIT("ModuleIO::read_dm","not ready to readin DM_R");
    }
#endif
    if(GlobalV::MY_RANK==0) ifs.close();

    GlobalV::ofs_running << " Finish reading density matrix." << std::endl;

    ModuleBase::timer::tick("ModuleIO","read_dm");
    return;
}
