#include "src_lcao/local_orbital_charge.h"
#include "src_pw/global.h"
#include "module_base/blas_connector.h"


void Local_Orbital_Charge::read_dm(const int &is, const string &fn)
{
    TITLE("Local_Orbital_Charge","read_dm");
    timer::tick("Local_Orbital_Charge","read_dm");

    ofs_running << "\n processor 0 is reading density matrix from file < " << fn << " > " << endl;
    //xiaohui modify 2015-03-25
    //bool quit_mesia = false;
    bool quit_abacus = false;

    ifstream ifs;
    if(MY_RANK==0)
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

            string name;
            ifs >> name;

            // check lattice constant, unit is Angstrom
            CHECK_DOUBLE(ifs,ucell.lat0 * BOHR_TO_A,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e11,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e12,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e13,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e21,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e22,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e23,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e31,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e32,quit);
            CHECK_DOUBLE(ifs,ucell.latvec.e33,quit);

            for(int it=0; it<ucell.ntype; it++)
            {
                CHECK_STRING(ifs,ucell.atoms[it].label,quit);
            }

            for(int it=0; it<ucell.ntype; it++)
            {
                CHECK_DOUBLE(ifs,ucell.atoms[it].na,quit);
            }

            string coordinate;
            ifs >> coordinate;

            for(int it=0; it<ucell.ntype; it++)
            {
                for(int ia=0; ia<ucell.atoms[it].na; ia++)
                {
                    CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].x,quit);
                    CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].y,quit);
                    CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].z,quit);
                }
            }

            CHECK_INT(ifs, NSPIN);
            if(NSPIN == 1||NSPIN == 4)
            {
                READ_VALUE(ifs, en.ef);
                ofs_running << " read in fermi energy = " << en.ef << endl;
            }
            else if(NSPIN == 2)
            {
                if(is==0)READ_VALUE(ifs, en.ef_up);
                else if(is==1)READ_VALUE(ifs, en.ef_dw);
            }
            else
            {
                WARNING_QUIT("read_dm","check nspin!");
            }
            CHECK_INT(ifs, NLOCAL);
            CHECK_INT(ifs, NLOCAL);
        }// If file exist, read in data.
    } // Finish reading the first part of density matrix.


#ifndef __MPI
    ofs_running << " Read SPIN = " << is+1 << " density matrix now." << endl;

    if(GAMMA_ONLY_LOCAL)
    {
        for(int i=0; i<NLOCAL; ++i)
        {
            for(int j=0; j<NLOCAL; ++j)
            {
                ifs >> DM[is][i][j];
            }
        }
    }
    else
    {
        WARNING_QUIT("Local_Orbital_Charge::read_dm","The nnrg should not be update");
        CHECK_INT(ifs,LNNR.nnrg);

        for(int i=0; i<LNNR.nnrg; ++i)
        {
            ifs >> DM_R[is][i];
        }
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
        WARNING_QUIT("Local_Orbital_Charge::read_dm","Can not find the density matrix file.");
    }


    if(NSPIN==1||NSPIN==4)
    {
        Parallel_Common::bcast_double(en.ef);
    }
    else if(NSPIN==2)
    {
        Parallel_Common::bcast_double(en.ef_up);
        Parallel_Common::bcast_double(en.ef_dw);
    }


    if(GAMMA_ONLY_LOCAL)
    {
        //ofs_running << " NLOCAL=" << NLOCAL << endl;
        //ofs_running << " lgd_now=" << lgd_now << endl;
        //ofs_running << " GridT.lgd=" << GridT.lgd << endl;

        double *tmp = new double[NLOCAL];
        for(int i=0; i<NLOCAL; ++i)
        {
            //ofs_running << " i=" << i << endl;
            ZEROS(tmp, NLOCAL);
            if(MY_RANK==0)
            {
                for(int j=0; j<NLOCAL; ++j)
                {
                    ifs >> tmp[j];
                }
            }
            Parallel_Common::bcast_double(tmp, NLOCAL);

            const int mu = GridT.trace_lo[i];
            if(mu >= 0)
            {   
                for(int j=0; j<NLOCAL; ++j)
                {
                    const int nu = GridT.trace_lo[j];
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
        WARNING_QUIT("Local_Orbital_Charge::read_dm","not ready to readin DM_R");
    }
#endif
    if(MY_RANK==0) ifs.close();

    ofs_running << " Finish reading density matrix." << endl;

    timer::tick("Local_Orbital_Charge","read_dm");
    return;
}
