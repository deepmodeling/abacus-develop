#include "FORCE_gamma.h"
#include "../src_pw/global.h"
#include "../module_base/timer.h"

void Force_LCAO_gamma::cal_fvl_dphi(
	ModuleBase::matrix& dm2d, 
	const bool isforce, 
    const bool isstress,
    ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& svl_dphi)
{   
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
	
	if(isforce)
	{
		ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_x, this->ParaV->nloc);
		ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_y, this->ParaV->nloc);
		ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_z, this->ParaV->nloc);
	}
    if(isstress)
    {
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_11, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_12, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_13, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_22, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_23, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_33, this->ParaV->nloc);
    }


    double* tmpDHx = new double[this->ParaV->nloc];
    double* tmpDHy = new double[this->ParaV->nloc];
    double* tmpDHz = new double[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS( tmpDHx, this->ParaV->nloc );
    ModuleBase::GlobalFunc::ZEROS( tmpDHy, this->ParaV->nloc );
    ModuleBase::GlobalFunc::ZEROS( tmpDHz, this->ParaV->nloc );
    for(int i=0; i<this->ParaV->nloc; ++i)
    {
        tmpDHx[i] = this->UHM->LM->DHloc_fixed_x[i];
        tmpDHy[i] = this->UHM->LM->DHloc_fixed_y[i];
        tmpDHz[i] = this->UHM->LM->DHloc_fixed_z[i];
        //std::cout << "  this->UHM->LM->DHloc_fixed_x=" <<  this->UHM->LM->DHloc_fixed_x[i] << std::endl;
        //std::cout << "  this->UHM->LM->DHloc_fixed_y=" <<  this->UHM->LM->DHloc_fixed_y[i] << std::endl;
        //std::cout << "  this->UHM->LM->DHloc_fixed_z=" <<  this->UHM->LM->DHloc_fixed_z[i] << std::endl;
    }

    //calculate dVL
    //calculate <dphi | VL | phi>

	// mohan add 2021, needs reconstruction!!!
    int istep = 1;
    GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_x, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_y, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_z, this->ParaV->nloc);

        for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        //  should not be set zero if VNA is used.
        //  ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixed_x,this->ParaV->nloc);
        //  ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixed_y,this->ParaV->nloc);
        //  ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixed_z,this->ParaV->nloc);
        this->UHM->GG.cal_force(GlobalC::pot.vr_eff1);


        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
            const int iat = GlobalC::ucell.iwt2iat[i];
            for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                //const int iat2 = GlobalC::ucell.iwt2iat[j];
                const int mu = this->ParaV->trace_loc_row[j];
                const int nu = this->ParaV->trace_loc_col[i];
                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * this->ParaV->ncol + nu;
                    //contribution from deriv of AO's in T+VNL term

                    double dm2d2 = 2.0 * dm2d(is, index);

					if(isforce)
					{
						fvl_dphi(iat,0) -= dm2d2 * ( this->UHM->LM->DHloc_fixed_x[index] + tmpDHx[index] );
						fvl_dphi(iat,1) -= dm2d2 * ( this->UHM->LM->DHloc_fixed_y[index] + tmpDHy[index] );
						fvl_dphi(iat,2) -= dm2d2 * ( this->UHM->LM->DHloc_fixed_z[index] + tmpDHz[index] );
					}
                    if(isstress)
                    {
                        svl_dphi(0,0) += dm2d2 * this->UHM->LM->DHloc_fixed_11[index];
                        svl_dphi(0,1) += dm2d2 * this->UHM->LM->DHloc_fixed_12[index];
                        svl_dphi(0,2) += dm2d2 * this->UHM->LM->DHloc_fixed_13[index];
                        svl_dphi(1,1) += dm2d2 * this->UHM->LM->DHloc_fixed_22[index];
                        svl_dphi(1,2) += dm2d2 * this->UHM->LM->DHloc_fixed_23[index];
                        svl_dphi(2,2) += dm2d2 * this->UHM->LM->DHloc_fixed_33[index];
                    }
                    //  std::cout << std::setw(5) << iat << std::setw(5) << iat2 
                    //  << std::setw(5) << mu << std::setw(5) << nu
                    //  << std::setw(15) << this->UHM->LM->DHloc_fixed_z[index] << std::endl;
                }
            }
        }

//          std::cout << "fvl_dphi:" << std::endl;
//          for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
//          {
//              std::cout << std::setw(5) << iat << std::setw(15) << fvl_dphi[iat][0] 
//              << std::setw(15) << fvl_dphi[iat][1]
//              << std::setw(15) << fvl_dphi[iat][2] << std::endl;
//          }


    } // end spin
    // test mohan tmp
//  test_gamma(this->UHM->LM->DHloc_fixed_x,"this->UHM->LM->DHloc_fixed_x");

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= GlobalC::ucell.omega;
            }
        }
    }

    delete[] tmpDHx;
    delete[] tmpDHy;
    delete[] tmpDHz;


    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
    return;
}



void Force_LCAO_gamma::cal_fvl_dphi(
	const std::vector<ModuleBase::matrix> &dm2d, 
	const bool isforce, 
    const bool isstress,
    ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& svl_dphi)
{   
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvl_dphi");

    ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_x, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_y, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_z, this->ParaV->nloc);
    if(GlobalV::CAL_STRESS)
    {
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_11, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_12, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_13, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_22, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_23, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_33, this->ParaV->nloc);
    }


    double* tmpDHx = new double[this->ParaV->nloc];
    double* tmpDHy = new double[this->ParaV->nloc];
    double* tmpDHz = new double[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS( tmpDHx, this->ParaV->nloc );
    ModuleBase::GlobalFunc::ZEROS( tmpDHy, this->ParaV->nloc );
    ModuleBase::GlobalFunc::ZEROS( tmpDHz, this->ParaV->nloc );
    for(int i=0; i<this->ParaV->nloc; ++i)
    {
        tmpDHx[i] = this->UHM->LM->DHloc_fixed_x[i];
        tmpDHy[i] = this->UHM->LM->DHloc_fixed_y[i];
        tmpDHz[i] = this->UHM->LM->DHloc_fixed_z[i];
        //std::cout << "  this->UHM->LM->DHloc_fixed_x=" <<  this->UHM->LM->DHloc_fixed_x[i] << std::endl;
        //std::cout << "  this->UHM->LM->DHloc_fixed_y=" <<  this->UHM->LM->DHloc_fixed_y[i] << std::endl;
        //std::cout << "  this->UHM->LM->DHloc_fixed_z=" <<  this->UHM->LM->DHloc_fixed_z[i] << std::endl;
    }

    //calculate dVL
    //calculate <dphi | VL | phi>

    int istep = 1;
    GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_x, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_y, this->ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS (this->UHM->LM->DHloc_fixed_z, this->ParaV->nloc);

        for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        //should not be set zero if VNA is used.
        //ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixed_x,this->ParaV->nloc);
        //ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixed_y,this->ParaV->nloc);
        //ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixed_z,this->ParaV->nloc);
        this->UHM->GG.cal_force(GlobalC::pot.vr_eff1);

        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
            const int iat = GlobalC::ucell.iwt2iat[i];
            for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                //const int iat2 = GlobalC::ucell.iwt2iat[j];
                const int mu = this->ParaV->trace_loc_row[j];
                const int nu = this->ParaV->trace_loc_col[i];
                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * this->ParaV->ncol + nu;
                    //contribution from deriv of AO's in T+VNL term

                    double dm2d2 = 2.0 * dm2d[is](nu, mu);

                    if(isforce)
					{
						fvl_dphi(iat,0) -= dm2d2 * ( this->UHM->LM->DHloc_fixed_x[index] + tmpDHx[index] );
						fvl_dphi(iat,1) -= dm2d2 * ( this->UHM->LM->DHloc_fixed_y[index] + tmpDHy[index] );
						fvl_dphi(iat,2) -= dm2d2 * ( this->UHM->LM->DHloc_fixed_z[index] + tmpDHz[index] );
					}
                    if(isstress)
                    {
                        svl_dphi(0,0) += dm2d2 * this->UHM->LM->DHloc_fixed_11[index];
                        svl_dphi(0,1) += dm2d2 * this->UHM->LM->DHloc_fixed_12[index];
                        svl_dphi(0,2) += dm2d2 * this->UHM->LM->DHloc_fixed_13[index];
                        svl_dphi(1,1) += dm2d2 * this->UHM->LM->DHloc_fixed_22[index];
                        svl_dphi(1,2) += dm2d2 * this->UHM->LM->DHloc_fixed_23[index];
                        svl_dphi(2,2) += dm2d2 * this->UHM->LM->DHloc_fixed_33[index];
                    }
                    //std::cout << std::setw(5) << iat << std::setw(5) << iat2 
                    //<< std::setw(5) << mu << std::setw(5) << nu
                    //<< std::setw(15) << this->UHM->LM->DHloc_fixed_z[index] << std::endl;
                }
            }
        }

        //std::cout << "fvl_dphi:" << std::endl;
        //for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
        //{
        //std::cout << std::setw(5) << iat << std::setw(15) << fvl_dphi[iat][0] 
        //<< std::setw(15) << fvl_dphi[iat][1]
        //<< std::setw(15) << fvl_dphi[iat][2] << std::endl;
        //}

    } // end spin
    // test mohan tmp
    //test_gamma(this->UHM->LM->DHloc_fixed_x,"this->UHM->LM->DHloc_fixed_x");

    if(isstress)
    {
        StressTools::stress_fill(1.0, GlobalC::ucell.omega, svl_dphi);
    }

    delete[] tmpDHx;
    delete[] tmpDHy;
    delete[] tmpDHz;

    ModuleBase::timer::tick("Force_LCAO_gamma", "cal_fvl_dphi");
    return;
}
