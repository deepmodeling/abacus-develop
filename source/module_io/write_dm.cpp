#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_io/dm_io.h"

//-------------------------------------------------
// NOTE for ModuleIO::write_dm
// I will give an example here, suppose we have a 4*4 
// density matrix (symmetry) which is
// 1.1 2.3 3.6 4.2
// 2.3 5.2 7.1 8.9
// 3.6 7.1 1.3 3.2  
// 4.2 8.9 3.2 2.4
// we use two processors, each one has 3 orbitals
// processor 1 has orbital index 1, 2, 4:
// ('no' means no value on this processor)
// 1.1 2.3 no  4.2
// 2.3 5.2 no  8.9
// no  no  no  no   
// 4.2 8.9 no  2.4
// processor 2 has orbital index 1, 3, 4;
// 1.1 no  3.6 4.2
// no  no  no  no 
// 3.6 no  1.3 3.2  
// 4.2 no  3.2 2.4
// now we want to reduce them and print out,
// we plan to reduce them one row by one row,
// then for the first row, we need to set the
// temparary array to 4 (GlobalV::NLOCAL in code),
// then we reduce first row, it will become
// 2.2 2.3 3.6 8.4,
// we notice that the first element and fourth
// element is doubled, that's because the density
// may have overlap, so we need to first count
// for each element, how many times it is duplicated
// on other processors, this is why there is
// a 'count' integer array in the code.
// UPDATED BY MOHAN 2014-05-18
void ModuleIO::write_dm(
#ifdef __MPI
	const int* trace_lo,
#endif
	const int &is, 
	const int &iter, 
	const std::string &fn, 
	int precision,
	int out_dm,
	double*** DM,
	const double& ef,
	const UnitCell* ucell,
    const int my_rank,
    const int nspin,
    const int nlocal)
{
    ModuleBase::TITLE("ModuleIO","write_dm");

	if (out_dm==0)
	{
		return;
	}

	ModuleBase::timer::tick("ModuleIO","write_dm");

	time_t start, end;
	std::ofstream ofs;

	if(my_rank==0)
	{
		start = time(NULL);

		ofs.open(fn.c_str());
		if (!ofs)
		{
			ModuleBase::WARNING("ModuleIO::write_dm","Can't create DENSITY MATRIX File!");
		}

		//GlobalV::ofs_running << "\n Output charge file." << std::endl;

		ofs << ucell->latName << std::endl;//1
		ofs << " " << ucell->lat0 * ModuleBase::BOHR_TO_A << std::endl;
		ofs << " " << ucell->latvec.e11 << " " << ucell->latvec.e12 << " " << ucell->latvec.e13 << std::endl;
		ofs << " " << ucell->latvec.e21 << " " << ucell->latvec.e22 << " " << ucell->latvec.e23 << std::endl;
		ofs << " " << ucell->latvec.e31 << " " << ucell->latvec.e32 << " " << ucell->latvec.e33 << std::endl;
		for(int it=0; it<ucell->ntype; it++)
		{
			ofs << " " << ucell->atoms[it].label;
		}
		ofs << std::endl;
		for(int it=0; it<ucell->ntype; it++)
		{
			ofs << " " << ucell->atoms[it].na;
		}
		ofs << std::endl;
		ofs << "Direct" << std::endl;

		for(int it=0; it<ucell->ntype; it++)
		{
			Atom* atom = &ucell->atoms[it];
			ofs << std::setprecision(15);
			for(int ia=0; ia<ucell->atoms[it].na; ia++)
			{
				ofs << " " << atom->taud[ia].x
					<< " " << atom->taud[ia].y
					<< " " << atom->taud[ia].z << std::endl;
			}
		}

		ofs << "\n " << nspin;
		ofs << "\n " << ef << " (fermi energy)";

		ofs << "\n  " << nlocal << " " << nlocal << std::endl;

		ofs << std::setprecision(precision);
		ofs << std::scientific;

	}

    //ofs << "\n " << GlobalV::GAMMA_ONLY_LOCAL << " (GAMMA ONLY LOCAL)" << std::endl;
#ifndef __MPI

    for(int i=0; i<nlocal; ++i)
    {
        for(int j=0; j<nlocal; ++j)
        {
            if(j%8==0) ofs << "\n";
            ofs << " " << DM[is][i][j];
        }
    }

#else
    //xiaohui modify 2014-06-18
    
    double* tmp = new double[nlocal];
    int* count = new int[nlocal];
    for (int i=0; i<nlocal; ++i)
    {
        // when reduce, there may be 'redundance', we need to count them.
        ModuleBase::GlobalFunc::ZEROS(count, nlocal);
        const int mu = trace_lo[i];
        if (mu >= 0)
        {
            for (int j=0; j<nlocal; ++j)
            {
                const int nu = trace_lo[j];
                if (nu >= 0)
                {
                    count[j]=1; 
                }
            }
        }
        Parallel_Reduce::reduce_all(count, nlocal);

        // reduce the density matrix for 'i' line.
        ModuleBase::GlobalFunc::ZEROS(tmp, nlocal);
        if (mu >= 0)
        {
            for (int j=0; j<nlocal; j++)
            {
                const int nu = trace_lo[j];
                if (nu >=0)
                {
                    tmp[j] = DM[is][mu][nu];
                    //GlobalV::ofs_running << " dmi=" << i << " j=" << j << " " << DM[is][mu][nu] << std::endl;
                }
            }
        }
        Parallel_Reduce::reduce_all(tmp, nlocal);

        if(my_rank==0)
        {
            for (int j=0; j<nlocal; j++)
            {
                if(j%8==0) ofs << "\n";
                if(count[j]>0)
                {
                    ofs << " " << tmp[j]/(double)count[j];
                }
                else
                {
                    ofs << " 0"; 
                }
            }
        }
    }
    delete[] tmp;
    delete[] count;
#endif
	if(my_rank==0)
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_dm",start,end);
		ofs.close();
	}
	ModuleBase::timer::tick("ModuleIO","write_dm");

    return;
}
