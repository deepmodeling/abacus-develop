#include "src_lcao/local_orbital_charge.h"
#include "src_pw/global.h"
#include "src_global/blas_connector.h"
#include "src_lcao/LCAO_nnr.h" // mohan add 2021-03-07

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
// temparary array to 4 (NLOCAL in code),
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
	const string &fn, 
	const int &precision)
{
    TITLE("Local_Orbital_Charge","write_dm");

	if (out_dm==0)
	{
		return;
	}
	else if(iter % out_dm != 0)
	{
		return; 
	}
	timer::tick("Local_Orbital_Charge","write_dm");

	time_t start, end;
	ofstream ofs;

	if(MY_RANK==0)
	{
		start = time(NULL);

		ofs.open(fn.c_str());
		if (!ofs)
		{
			WARNING("Charge::write_rho","Can't create Charge File!");
		}

		//ofs_running << "\n Output charge file." << endl;

		ofs << ucell.latName << endl;//1
		ofs << " " << ucell.lat0 * BOHR_TO_A << endl;
		ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;
		for(int it=0; it<ucell.ntype; it++)
		{
			ofs << " " << ucell.atoms[it].label;
		}
		ofs << endl;
		for(int it=0; it<ucell.ntype; it++)
		{
			ofs << " " << ucell.atoms[it].na;
		}
		ofs << endl;
		ofs << "Direct" << endl;

		for(int it=0; it<ucell.ntype; it++)
		{
			Atom* atom = &ucell.atoms[it];
			ofs << setprecision(15);
			for(int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				ofs << " " << atom->taud[ia].x
					<< " " << atom->taud[ia].y
					<< " " << atom->taud[ia].z << endl;
			}
		}

		ofs << "\n " << NSPIN;
		if(NSPIN==1||NSPIN==4)
		{
			ofs << "\n " << en.ef << " (fermi energy)";
		}
		else if(NSPIN==2)
		{
			if(is==0)ofs << "\n " << en.ef_up << " (fermi energy for spin=1)";
			else if(is==1)ofs << "\n " << en.ef_dw << " (fermi energy for spin=2)";
		}
		else
		{
			WARNING_QUIT("write_rho","check nspin!");
		}


		ofs << "\n  " << NLOCAL << " " << NLOCAL << endl;

		ofs << setprecision(precision);
		ofs << scientific;

	}

    //ofs << "\n " << GAMMA_ONLY_LOCAL << " (GAMMA ONLY LOCAL)" << endl;
#ifndef __MPI
    if(GAMMA_ONLY_LOCAL)
    {
        for(int i=0; i<NLOCAL; ++i)
        {
            for(int j=0; j<NLOCAL; ++j)
            {
                if(j%8==0) ofs << "\n";
                ofs << " " << this->DM[is][i][j];
            }
        }
    }
    else
    {
        WARNING_QUIT("write_dm","not ready yet");
        ofs << " " << LNNR.nnrg << " (nnrg)" << endl;
        for(int i=0; i<LNNR.nnrg; ++i)
        {
            if(i%8==0) ofs << "\n";
            ofs << " " << this->DM_R[is][i];
        }
    }
#else
    if(GAMMA_ONLY_LOCAL)
    {
        //xiaohui modify 2014-06-18
        
        double* tmp = new double[NLOCAL];
        int* count = new int[NLOCAL];
        for (int i=0; i<NLOCAL; ++i)
        {
            // when reduce, there may be 'redundance', we need to count them.
            ZEROS(count, NLOCAL);
            const int mu = GridT.trace_lo[i];
            if (mu >= 0)
            {
                for (int j=0; j<NLOCAL; ++j)
                {
                    const int nu = GridT.trace_lo[j];
                    if (nu >= 0)
                    {
                        count[j]=1; 
                    }
                }
            }
            Parallel_Reduce::reduce_int_all( count, NLOCAL );

            // reduce the density matrix for 'i' line.
            ZEROS(tmp, NLOCAL);
            if (mu >= 0)
            {
                for (int j=0; j<NLOCAL; j++)
                {
                    const int nu = GridT.trace_lo[j];
                    if (nu >=0)
                    {
                        tmp[j] = DM[is][mu][nu];
                        //ofs_running << " dmi=" << i << " j=" << j << " " << DM[is][mu][nu] << endl;
                    }
                }
            }
            Parallel_Reduce::reduce_double_all( tmp, NLOCAL );

            if(MY_RANK==0)
            {
                for (int j=0; j<NLOCAL; j++)
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
        //for(int i=0; i<NLOCAL; ++i)
        //{
        //  for(int j=0; j<NLOCAL; ++j)
        //  {
        //      if(j%8==0) ofs << "\n";
        //      ofs << " " << this->DM[is][i][j];
        //  }
        //}

    }
    else
    {
        ofs << " " << LNNR.nnrg << " (nnrg)" << endl;
        WARNING_QUIT("local_orbital_charge","not ready to output DM_R");
    }
#endif
	if(MY_RANK==0)
	{
		end = time(NULL);
		OUT_TIME("write_rho",start,end);
		ofs.close();
	}
	timer::tick("Local_Orbital_Charge","write_dm");

    return;
}
