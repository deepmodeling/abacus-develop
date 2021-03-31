#include "FORCE_gamma.h"
#include "../src_pw/global.h"

void Force_LCAO_gamma::cal_fvl_dphi(
	matrix& dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& fvl_dphi, 
	matrix& svl_dphi)
{   
    TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    timer::tick("Force_LCAO_gamma","cal_fvl_dphi",'G');
	
	if(isforce)
	{
		ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
		ZEROS (LM.DHloc_fixed_z, ParaO.nloc);
	}
    if(isstress)
    {
        ZEROS (LM.DHloc_fixed_11, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_12, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_13, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_22, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_23, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_33, ParaO.nloc);
    }


    double* tmpDHx = new double[ParaO.nloc];
    double* tmpDHy = new double[ParaO.nloc];
    double* tmpDHz = new double[ParaO.nloc];
    ZEROS( tmpDHx, ParaO.nloc );
    ZEROS( tmpDHy, ParaO.nloc );
    ZEROS( tmpDHz, ParaO.nloc );
    for(int i=0; i<ParaO.nloc; ++i)
    {
        tmpDHx[i] = LM.DHloc_fixed_x[i];
        tmpDHy[i] = LM.DHloc_fixed_y[i];
        tmpDHz[i] = LM.DHloc_fixed_z[i];
        //cout << "  LM.DHloc_fixed_x=" <<  LM.DHloc_fixed_x[i] << endl;
        //cout << "  LM.DHloc_fixed_y=" <<  LM.DHloc_fixed_y[i] << endl;
        //cout << "  LM.DHloc_fixed_z=" <<  LM.DHloc_fixed_z[i] << endl;
    }

    //calculate dVL
    //calculate <dphi | VL | phi>

	// mohan add 2021, needs reconstruction!!!
    int istep = 1;
    pot.init_pot(istep, pw.strucFac);

    for(int is=0; is<NSPIN; ++is)
    {
        CURRENT_SPIN = is;

        ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_z, ParaO.nloc);

        for(int ir=0; ir<pw.nrxx; ++ir)
        {
            pot.vr_eff1[ir] = pot.vr_eff(CURRENT_SPIN, ir);
        }

        //  should not be set zero if VNA is used.
        //  ZEROS(LM.DHloc_fixed_x,ParaO.nloc);
        //  ZEROS(LM.DHloc_fixed_y,ParaO.nloc);
        //  ZEROS(LM.DHloc_fixed_z,ParaO.nloc);
        UHM.GG.cal_force(pot.vr_eff1);


        for(int i=0; i<NLOCAL; i++)
        {
            const int iat = ucell.iwt2iat[i];
            for(int j=0; j<NLOCAL; j++)
            {
                const int iat2 = ucell.iwt2iat[j];
                const int mu = ParaO.trace_loc_row[j];
                const int nu = ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * ParaO.ncol + nu;
                    //contribution from deriv of AO's in T+VNL term

                    double dm2d2 = 2.0 * dm2d(is, index);

					if(isforce)
					{
						fvl_dphi(iat,0) -= dm2d2 * ( LM.DHloc_fixed_x[index] + tmpDHx[index] );
						fvl_dphi(iat,1) -= dm2d2 * ( LM.DHloc_fixed_y[index] + tmpDHy[index] );
						fvl_dphi(iat,2) -= dm2d2 * ( LM.DHloc_fixed_z[index] + tmpDHz[index] );
					}
                    if(isstress)
                    {
                        svl_dphi(0,0) += dm2d2 * LM.DHloc_fixed_11[index];
                        svl_dphi(0,1) += dm2d2 * LM.DHloc_fixed_12[index];
                        svl_dphi(0,2) += dm2d2 * LM.DHloc_fixed_13[index];
                        svl_dphi(1,1) += dm2d2 * LM.DHloc_fixed_22[index];
                        svl_dphi(1,2) += dm2d2 * LM.DHloc_fixed_23[index];
                        svl_dphi(2,2) += dm2d2 * LM.DHloc_fixed_33[index];
                    }
                    //  cout << setw(5) << iat << setw(5) << iat2 
                    //  << setw(5) << mu << setw(5) << nu
                    //  << setw(15) << LM.DHloc_fixed_z[index] << endl;
                }
            }
        }

//          cout << "fvl_dphi:" << endl;
//          for(int iat=0; iat<ucell.nat; ++iat)
//          {
//              cout << setw(5) << iat << setw(15) << fvl_dphi[iat][0] 
//              << setw(15) << fvl_dphi[iat][1]
//              << setw(15) << fvl_dphi[iat][2] << endl;
//          }


    } // end spin
    // test mohan tmp
//  test_gamma(LM.DHloc_fixed_x,"LM.DHloc_fixed_x");

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= ucell.omega;
            }
        }
    }

    delete[] tmpDHx;
    delete[] tmpDHy;
    delete[] tmpDHz;


    timer::tick("Force_LCAO_gamma","cal_fvl_dphi",'G');
    return;
}



void Force_LCAO_gamma::cal_fvl_dphi(
	const std::vector<matrix> &dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& fvl_dphi, 
	matrix& svl_dphi)
{   
    TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    timer::tick("Force_LCAO_gamma","cal_fvl_dphi",'G');

    ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_z, ParaO.nloc);
    if(STRESS)
    {
        ZEROS (LM.DHloc_fixed_11, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_12, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_13, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_22, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_23, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_33, ParaO.nloc);
    }


    double* tmpDHx = new double[ParaO.nloc];
    double* tmpDHy = new double[ParaO.nloc];
    double* tmpDHz = new double[ParaO.nloc];
    ZEROS( tmpDHx, ParaO.nloc );
    ZEROS( tmpDHy, ParaO.nloc );
    ZEROS( tmpDHz, ParaO.nloc );
    for(int i=0; i<ParaO.nloc; ++i)
    {
        tmpDHx[i] = LM.DHloc_fixed_x[i];
        tmpDHy[i] = LM.DHloc_fixed_y[i];
        tmpDHz[i] = LM.DHloc_fixed_z[i];
        //cout << "  LM.DHloc_fixed_x=" <<  LM.DHloc_fixed_x[i] << endl;
        //cout << "  LM.DHloc_fixed_y=" <<  LM.DHloc_fixed_y[i] << endl;
        //cout << "  LM.DHloc_fixed_z=" <<  LM.DHloc_fixed_z[i] << endl;
    }

    //calculate dVL
    //calculate <dphi | VL | phi>

    int istep = 1;
    pot.init_pot(istep, pw.strucFac);

    for(int is=0; is<NSPIN; ++is)
    {
        CURRENT_SPIN = is;

        ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_z, ParaO.nloc);

        for(int ir=0; ir<pw.nrxx; ++ir)
        {
            pot.vr_eff1[ir] = pot.vr_eff(CURRENT_SPIN, ir);
        }

        //should not be set zero if VNA is used.
        //ZEROS(LM.DHloc_fixed_x,ParaO.nloc);
        //ZEROS(LM.DHloc_fixed_y,ParaO.nloc);
        //ZEROS(LM.DHloc_fixed_z,ParaO.nloc);
        UHM.GG.cal_force(pot.vr_eff1);

        for(int i=0; i<NLOCAL; i++)
        {
            const int iat = ucell.iwt2iat[i];
            for(int j=0; j<NLOCAL; j++)
            {
                const int iat2 = ucell.iwt2iat[j];
                const int mu = ParaO.trace_loc_row[j];
                const int nu = ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * ParaO.ncol + nu;
                    //contribution from deriv of AO's in T+VNL term

                    double dm2d2 = 2.0 * dm2d[is](nu, mu);

                    if(isforce)
					{
						fvl_dphi(iat,0) -= dm2d2 * ( LM.DHloc_fixed_x[index] + tmpDHx[index] );
						fvl_dphi(iat,1) -= dm2d2 * ( LM.DHloc_fixed_y[index] + tmpDHy[index] );
						fvl_dphi(iat,2) -= dm2d2 * ( LM.DHloc_fixed_z[index] + tmpDHz[index] );
					}
                    if(isstress)
                    {
                        svl_dphi(0,0) += dm2d2 * LM.DHloc_fixed_11[index];
                        svl_dphi(0,1) += dm2d2 * LM.DHloc_fixed_12[index];
                        svl_dphi(0,2) += dm2d2 * LM.DHloc_fixed_13[index];
                        svl_dphi(1,1) += dm2d2 * LM.DHloc_fixed_22[index];
                        svl_dphi(1,2) += dm2d2 * LM.DHloc_fixed_23[index];
                        svl_dphi(2,2) += dm2d2 * LM.DHloc_fixed_33[index];
                    }
                    //cout << setw(5) << iat << setw(5) << iat2 
                    //<< setw(5) << mu << setw(5) << nu
                    //<< setw(15) << LM.DHloc_fixed_z[index] << endl;
                }
            }
        }

        //cout << "fvl_dphi:" << endl;
        //for(int iat=0; iat<ucell.nat; ++iat)
        //{
        //cout << setw(5) << iat << setw(15) << fvl_dphi[iat][0] 
        //<< setw(15) << fvl_dphi[iat][1]
        //<< setw(15) << fvl_dphi[iat][2] << endl;
        //}

    } // end spin
    // test mohan tmp
    //test_gamma(LM.DHloc_fixed_x,"LM.DHloc_fixed_x");

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= ucell.omega;
            }
        }
    }

    delete[] tmpDHx;
    delete[] tmpDHy;
    delete[] tmpDHz;

    timer::tick("Force_LCAO_gamma","cal_fvl_dphi",'G');
    return;
}
