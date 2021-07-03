#include "global.h"
#include "threshold_elec.h"

Threshold_Elec::Threshold_Elec() 
{ 
	dr2 = 0.0; 
	conv_elec = false; 
}

void Threshold_Elec::set_ethr(void) const
{
	TITLE("Threshold_Elec","set_ethr");
    //========================================================================
    // setup ethr see setup.f90 Page 5/14
    // setup ethr, the convergence threshold for eigenvalues
    // tr2: the convergence threshold for self-consistency (rho or potential);
    //========================================================================

    //===================
    // nscf calculations
    //===================
    if (CALCULATION=="nscf")
    {
        if (ETHR == 0.0)
        {
            ETHR = 0.1 * std::min(1.0e-2, DRHO2 / CHR.nelec);
        }
    }
    //=================
    // self consistent
    //=================
    else if(CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax" || CALCULATION=="scf-sto")//qianrui 2021-2-20
    {
        if (ETHR == 0.0)
        {
            if (pot.start_pot == "file")
            {
                //======================================================
                // if you think that the starting potential is good
                // do not spoil it with a louly first diagonalization:
                // set a strict ethr in the input file ()diago_the_init
                //======================================================
                ETHR = 1.0e-5;
            }
            else
            {
                //=======================================================
                // starting atomic potential is probably far from scf
                // don't waste iterations in the first diagonalization
                //=======================================================
                ETHR = 1.0e-2;
            }
        }
    }

    return;
}

void Threshold_Elec::update_ethr(const int &iter)
{
    //======================================================
    // Convergence threshold for iterative diagonalization
    // is automatically updated during self consistency
    //======================================================
    if (iter > 1)
    {
        if (iter == 2)
        {
            ETHR = 1.e-2;
        }

		//----------------------------
		// ETHR changes in CG Method. 
		// mohan update 2012-03-26
		// mohan update 2012-02-08
		//----------------------------
		if(BASIS_TYPE=="lcao")
		{
			ETHR = std::min( ETHR, 0.01*dr2/ std::max(1.0, CHR.nelec));
		}
		// mohan update 2009-09-04
		else
		{
			ETHR = std::min( ETHR, 0.1*dr2/ std::max(1.0, CHR.nelec));
			//cout << " new ethr = " << ETHR << endl;
		}

    }
    return;
}

void Threshold_Elec::iter_end(ofstream &ofs)
{
	if(OUT_LEVEL != "m") 
	{
		print_eigenvalue(ofs);
	}
    return;
}


void Threshold_Elec::print_eigenvalue(ofstream &ofs)
{
	bool wrong = false;
	for(int ik=0; ik<kv.nks; ++ik)
	{
		for(int ib=0; ib<NBANDS; ++ib)
		{
			if( abs( wf.ekb[ik][ib] ) > 1.0e10)
			{
				ofs_warning << " ik=" << ik+1 << " ib=" << ib+1 << " " << wf.ekb[ik][ib] << " Ry" << endl;
				wrong = true;
			}
		}
	}
	if(wrong)
	{
		WARNING_QUIT("Threshold_Elec::print_eigenvalue","Eigenvalues are too large!");
	}


	if(MY_RANK!=0) 
	{
		return;
	}

	TITLE("Threshold_Elec","print_eigenvalue");

    ofs << "\n STATE ENERGY(eV) AND OCCUPATIONS.";
	ofs << setprecision(5);
    for (int ik = 0;ik < kv.nks;ik++)
    {
        if (NSPIN==2)
        {
            if (kv.isk[ik] == 0)ofs << "\n spin up :";
            if (kv.isk[ik] == 1)ofs << "\n spin down :";
        }
        
        if (NSPIN==2)
        {
            if (kv.isk[ik] == 0)
            {
                ofs << " " << ik+1 << "/" << kv.nks/2 << " kpoint (Cartesian) = "
                << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z
                << " (" << kv.ngk[ik] << " pws)" << endl;

                ofs << setprecision(6);

            }
            if (kv.isk[ik] == 1)
            {
                ofs << " " << ik+1-kv.nks/2 << "/" << kv.nks/2 << " kpoint (Cartesian) = "
                << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z
                << " (" << kv.ngk[ik] << " pws)" << endl;

                ofs << setprecision(6);

			}
		}       // Pengfei Li  added  14-9-9
		else	
		{
			ofs << " " << ik+1 << "/" << kv.nks << " kpoint (Cartesian) = " 
				<< kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z 
				<< " (" << kv.ngk[ik] << " pws)" << endl; 

			ofs << setprecision(6);
		}

		//----------------------
		// no energy to output
		//----------------------
		if(KS_SOLVER=="selinv")
		{
			ofs << " USING SELINV, NO BAND ENERGY IS AVAILABLE." << endl;
		}
		//----------------------
		// output energy
		//----------------------
		else
		{
			//ofs << setw(12) << kv.ngk[ik] << " PWs ";
			ofs_running << setprecision(6);
			ofs_running << setiosflags(ios::showpoint);
			for (int ib = 0; ib < NBANDS; ib++)
			{
				ofs << " [spin" << kv.isk[ik]+1 << "_state] " << setw(8) << ib+1 
				<< setw(15) << wf.ekb[ik][ib] * Ry_to_eV << setw(15) << wf.wg(ik, ib) << endl;
			}
			ofs << endl;
		}
    }//end ik
    return;
}
