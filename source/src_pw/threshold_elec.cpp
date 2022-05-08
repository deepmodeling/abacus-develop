#include "global.h"
#include "threshold_elec.h"

Threshold_Elec::Threshold_Elec()
{
	scf_thr = 0.0;
	conv_elec = false;
}

void Threshold_Elec::set_pw_diag_thr(void) const
{
	ModuleBase::TITLE("Threshold_Elec","set_pw_diag_thr");
    //========================================================================
    // setup pw_diag_thr see setup.f90 Page 5/14
    // setup pw_diag_thr, the convergence threshold for eigenvalues
    // tr2: the convergence threshold for self-consistency (rho or potential);
    //========================================================================

    //===================
    // nscf calculations
    //===================
    if (GlobalV::CALCULATION=="nscf")
    {
        if (abs(GlobalV::PW_DIAG_THR - 1.0e-2) < 1.0e-10)
        {
            GlobalV::PW_DIAG_THR = 0.1 * std::min(1.0e-2, GlobalV::SCF_THR / GlobalC::CHR.nelec);
        }
    }
    //=================
    // self consistent
    //=================
    else if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="scf-sto")//qianrui 2021-2-20
    {
        if (abs(GlobalV::PW_DIAG_THR - 1.0e-2) < 1.0e-10)
        {
            if (GlobalC::pot.init_chg == "file")
            {
                //======================================================
                // if you think that the starting potential is good
                // do not spoil it with a louly first diagonalization:
                // set a strict pw_diag_thr in the input file ()diago_the_init
                //======================================================
                GlobalV::PW_DIAG_THR = 1.0e-5;
            }
            else
            {
                //=======================================================
                // starting atomic potential is probably far from scf
                // don't waste iterations in the first diagonalization
                //=======================================================
                GlobalV::PW_DIAG_THR = 1.0e-2;
            }
        }
    }

    return;
}

void Threshold_Elec::update_pw_diag_thr(const int &iter)
{
    //======================================================
    // Convergence threshold for iterative diagonalization
    // is automatically updated during self consistency
    //======================================================
    if (iter > 1)
    {
        if (iter == 2)
        {
            GlobalV::PW_DIAG_THR = 1.e-2;
        }

		//----------------------------
		// GlobalV::PW_DIAG_THR changes in CG Method.
		// mohan update 2012-03-26
		// mohan update 2012-02-08
		//----------------------------
		if(GlobalV::BASIS_TYPE=="lcao")
		{
			GlobalV::PW_DIAG_THR = std::min( GlobalV::PW_DIAG_THR, 0.01*scf_thr/ std::max(1.0, GlobalC::CHR.nelec));
		}
		// mohan update 2009-09-04
		else
		{
			GlobalV::PW_DIAG_THR = std::min( GlobalV::PW_DIAG_THR, 0.1*scf_thr/ std::max(1.0, GlobalC::CHR.nelec));
			//std::cout << " new pw_diag_thr = " << GlobalV::PW_DIAG_THR << std::endl;
		}

    }
    else
    {
        if(GlobalV::CALCULATION=="md"||GlobalV::CALCULATION=="relax"||GlobalV::CALCULATION=="cell-relax")
        {
            GlobalV::PW_DIAG_THR = std::max(GlobalV::PW_DIAG_THR, INPUT.pw_diag_thr);
        }
    }
    return;
}

void Threshold_Elec::print_eigenvalue(std::ofstream &ofs)
{
	bool wrong = false;
	for(int ik=0; ik<GlobalC::kv.nks; ++ik)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ++ib)
		{
			if( abs( GlobalC::wf.ekb[ik][ib] ) > 1.0e10)
			{
				GlobalV::ofs_warning << " ik=" << ik+1 << " ib=" << ib+1 << " " << GlobalC::wf.ekb[ik][ib] << " Ry" << std::endl;
				wrong = true;
			}
		}
	}
	if(wrong)
	{
		ModuleBase::WARNING_QUIT("Threshold_Elec::print_eigenvalue","Eigenvalues are too large!");
	}


	if(GlobalV::MY_RANK!=0)
	{
		return;
	}

	ModuleBase::TITLE("Threshold_Elec","print_eigenvalue");

    ofs << "\n STATE ENERGY(eV) AND OCCUPATIONS ";
	ofs << std::setprecision(5);
    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {
        if(ik==0)
        {
            ofs << "   NSPIN == " << GlobalV::NSPIN << std::endl;
            if(GlobalV::NSPIN == 2)
            {
                ofs << "SPIN UP : " << std::endl;
            }
        }
        else if(ik == GlobalC::kv.nks/2)
        {
            if(GlobalV::NSPIN == 2)
            {
                ofs << "SPIN DOWN : " << std::endl;
            }
        }

        if (GlobalV::NSPIN==2)
        {
            if (GlobalC::kv.isk[ik] == 0)
            {
                ofs << " " << ik+1 << "/" << GlobalC::kv.nks/2 << " kpoint (Cartesian) = "
                << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);

            }
            if (GlobalC::kv.isk[ik] == 1)
            {
                ofs << " " << ik+1-GlobalC::kv.nks/2 << "/" << GlobalC::kv.nks/2 << " kpoint (Cartesian) = "
                << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);

			}
		}       // Pengfei Li  added  14-9-9
		else
		{
			ofs << " " << ik+1 << "/" << GlobalC::kv.nks << " kpoint (Cartesian) = "
				<< GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
				<< " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

			ofs << std::setprecision(6);
		}

		//----------------------
		// no energy to output
		//----------------------
		if(GlobalV::KS_SOLVER=="selinv")
		{
			ofs << " USING SELINV, NO BAND ENERGY IS AVAILABLE." << std::endl;
		}
		//----------------------
		// output energy
		//----------------------
		else
		{
			GlobalV::ofs_running << std::setprecision(6);
			GlobalV::ofs_running << std::setiosflags(ios::showpoint);
			for (int ib = 0; ib < GlobalV::NBANDS; ib++)
			{
				ofs << std::setw(8) << ib+1 
				    << std::setw(15) << GlobalC::wf.ekb[ik][ib] * ModuleBase::Ry_to_eV 
                    << std::setw(15) << GlobalC::wf.wg(ik, ib) << std::endl;
			}
			ofs << std::endl;
		}
    }//end ik
    return;
}
