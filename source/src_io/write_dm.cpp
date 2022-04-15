#include "../src_lcao/local_orbital_charge.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"

//-------------------------------------------------
// NOTE for Local_Orbital_Charge::write_dm
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
void Local_Orbital_Charge::write_dm(
	const int &is, 
	const int &iter, 
	const std::string &fn, 
	const int &precision)
{
    ModuleBase::TITLE("Local_Orbital_Charge","write_dm");

	if (out_dm==0)
	{
		return;
	}
	else if(iter % out_dm != 0)
	{
		return; 
	}
	ModuleBase::timer::tick("Local_Orbital_Charge","write_dm");

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

		//GlobalV::ofs_running << "\n Output charge file." << std::endl;

		ofs << GlobalC::ucell.latName << std::endl;//1
		ofs << " " << GlobalC::ucell.lat0 * ModuleBase::BOHR_TO_A << std::endl;
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
			Atom* atom = &GlobalC::ucell.atoms[it];
			ofs << std::setprecision(15);
			for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				ofs << " " << atom->taud[ia].x
					<< " " << atom->taud[ia].y
					<< " " << atom->taud[ia].z << std::endl;
			}
		}

		ofs << "\n " << GlobalV::NSPIN;
		if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
		{
			ofs << "\n " << GlobalC::en.ef << " (fermi energy)";
		}
		else if(GlobalV::NSPIN==2)
		{
			if(is==0)ofs << "\n " << GlobalC::en.ef_up << " (fermi energy for spin=1)";
			else if(is==1)ofs << "\n " << GlobalC::en.ef_dw << " (fermi energy for spin=2)";
		}
		else
		{
			ModuleBase::WARNING_QUIT("write_rho","check nspin!");
		}


		ofs << "\n  " << GlobalV::NLOCAL << " " << GlobalV::NLOCAL << std::endl;

		ofs << std::setprecision(precision);
		ofs << scientific;

	}

    //ofs << "\n " << GlobalV::GAMMA_ONLY_LOCAL << " (GAMMA ONLY LOCAL)" << std::endl;
#ifndef __MPI
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            for(int j=0; j<GlobalV::NLOCAL; ++j)
            {
                if(j%8==0) ofs << "\n";
                ofs << " " << this->DM[is][i][j];
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("write_dm","not ready yet");
        ofs << " " << GlobalC::GridT.nnrg << " (nnrg)" << std::endl;
        for(int i=0; i<GlobalC::GridT.nnrg; ++i)
        {
            if(i%8==0) ofs << "\n";
            ofs << " " << this->DM_R[is][i];
        }
    }
#else
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        //xiaohui modify 2014-06-18
        
        double* tmp = new double[GlobalV::NLOCAL];
        int* count = new int[GlobalV::NLOCAL];
        for (int i=0; i<GlobalV::NLOCAL; ++i)
        {
            // when reduce, there may be 'redundance', we need to count them.
            ModuleBase::GlobalFunc::ZEROS(count, GlobalV::NLOCAL);
            const int mu = GlobalC::GridT.trace_lo[i];
            if (mu >= 0)
            {
                for (int j=0; j<GlobalV::NLOCAL; ++j)
                {
                    const int nu = GlobalC::GridT.trace_lo[j];
                    if (nu >= 0)
                    {
                        count[j]=1; 
                    }
                }
            }
            Parallel_Reduce::reduce_int_all( count, GlobalV::NLOCAL );

            // reduce the density matrix for 'i' line.
            ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
            if (mu >= 0)
            {
                for (int j=0; j<GlobalV::NLOCAL; j++)
                {
                    const int nu = GlobalC::GridT.trace_lo[j];
                    if (nu >=0)
                    {
                        tmp[j] = DM[is][mu][nu];
                        //GlobalV::ofs_running << " dmi=" << i << " j=" << j << " " << DM[is][mu][nu] << std::endl;
                    }
                }
            }
            Parallel_Reduce::reduce_double_all( tmp, GlobalV::NLOCAL );

            if(GlobalV::MY_RANK==0)
            {
                for (int j=0; j<GlobalV::NLOCAL; j++)
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
        
        //xiaohui add 2014-06-18
        //for(int i=0; i<GlobalV::NLOCAL; ++i)
        //{
        //  for(int j=0; j<GlobalV::NLOCAL; ++j)
        //  {
        //      if(j%8==0) ofs << "\n";
        //      ofs << " " << this->DM[is][i][j];
        //  }
        //}

    }
    else
    {
        ofs << " " << GlobalC::GridT.nnrg << " (nnrg)" << std::endl;
        ModuleBase::WARNING_QUIT("local_orbital_charge","not ready to output DM_R");
    }
#endif
	if(GlobalV::MY_RANK==0)
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_rho",start,end);
		ofs.close();
	}
	ModuleBase::timer::tick("Local_Orbital_Charge","write_dm");

    return;
}
