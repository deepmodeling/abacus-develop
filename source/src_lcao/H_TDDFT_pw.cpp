#include "../src_pw/potential.h"
#include "../src_pw/global.h"
#include "ELEC_evolve.h"
#include "../module_base/timer.h"

//==========================================================
// this function aims to add external time-dependent potential
// (eg: linear potential) used in tddft
// fuxiang add in 2017-05
//==========================================================
void Potential::set_vrs_tddft(const int istep)
{
    ModuleBase::TITLE("Potential","set_vrs_tddft");
    ModuleBase::timer::tick("Potential","set_vrs_tddft");

    for (int is = 0;is < GlobalV::NSPIN;is++)
    {
        //====================================================
        // add external linear potential, fuxiang add in 2017/05
        //====================================================

        //const int timescale = 1;  // get the time that vext influences;
        if (istep >= ELEC_evolve::td_timescale)
        {
            for (int i = 0;i < GlobalC::pw.nrxx;i++)
            {
                this->vr_eff(is, i) = this->vltot[i] + this->vr(is, i);
            }
            std::cout << "vext = 0! " << std::endl;
        }
        else
        {
            this->vextold = new double[GlobalC::pw.nrxx];
            this->vext = new double[GlobalC::pw.nrxx];
            const int yz = GlobalC::pw.ncy*GlobalC::pw.nczp;
            int index, i, j, k;

            for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
            {
                index = ir;
                i     = index / yz; // get the z, z is the fastest
                index = index - yz * i;// get (x,y)
                j     = index / GlobalC::pw.nczp;// get y
                k     = index - GlobalC::pw.nczp*j + GlobalC::pw.nczp_start;// get x

                if(ELEC_evolve::td_vext_dire == 1)
                {
                    if (k<GlobalC::pw.ncx*0.05)
					{
						this->vextold[ir] = (0.019447*k/GlobalC::pw.ncx-0.001069585)*GlobalC::ucell.lat0;
					}
                    else if (k>=GlobalC::pw.ncx*0.05 && k<GlobalC::pw.ncx*0.95)
					{
						this->vextold[ir] = -0.0019447*k/GlobalC::pw.ncx*GlobalC::ucell.lat0;
					}
                    else if (k>=GlobalC::pw.ncx*0.95)
					{
						this->vextold[ir] = (0.019447*(1.0*k/GlobalC::pw.ncx-1)-0.001069585)*GlobalC::ucell.lat0;
					}
                }
                else if(ELEC_evolve::td_vext_dire == 2)
                {
                    if (j<GlobalC::pw.ncx*0.05)
					{
						this->vextold[ir] = (0.019447*j/GlobalC::pw.ncx-0.001069585)*GlobalC::ucell.lat0;
					}
                    else if (j>=GlobalC::pw.ncx*0.05 && j<GlobalC::pw.ncx*0.95)
					{
						this->vextold[ir] = -0.0019447*j/GlobalC::pw.ncx*GlobalC::ucell.lat0;
					}
                    else if (j>=GlobalC::pw.ncx*0.95)
					{
						this->vextold[ir] = (0.019447*(1.0*j/GlobalC::pw.ncx-1)-0.001069585)*GlobalC::ucell.lat0;
					}
                }
                else if(ELEC_evolve::td_vext_dire == 3)
                {
                    if (i<GlobalC::pw.ncx*0.05)
					{
						this->vextold[ir] = (0.019447*i/GlobalC::pw.ncx-0.001069585)*GlobalC::ucell.lat0;
					}
                    else if (i>=GlobalC::pw.ncx*0.05 && i<GlobalC::pw.ncx*0.95)
					{
						this->vextold[ir] = -0.0019447*i/GlobalC::pw.ncx*GlobalC::ucell.lat0;
					}
                    else if (i>=GlobalC::pw.ncx*0.95)
					{
						this->vextold[ir] = (0.019447*(1.0*i/GlobalC::pw.ncx-1)-0.001069585)*GlobalC::ucell.lat0;
					}
                }

                // Gauss
				if (ELEC_evolve::td_vexttype == 1)
				{
					const double w = 22.13;    // eV
					const double sigmasquare = 700;
					const double timecenter = 700;
					//Notice: these three parameters should be written in INPUT. I will correct soon.
					const double timenow = (istep-timecenter)*ELEC_evolve::td_scf_thr*41.34; //41.34 is conversion factor of fs-a.u.
					this->vext[ir] = this->vextold[ir]*cos(w/27.2116*timenow)*exp(-timenow*timenow*0.5/(sigmasquare))*0.25;  //0.1 is modified in 2018/1/12
				}

                //HHG of H atom
				if (ELEC_evolve::td_vexttype == 2)
				{
					const double w_h = 0.0588; //a.u.
					const int stepcut1 = 1875;
					const int stepcut2 = 5625;
					const int stepcut3 = 7500;
					// The parameters should be written in INPUT!
					if(istep < stepcut1)
					{
						this->vext[ir] = this->vextold[ir]*2.74*istep/stepcut1*cos(w_h*istep*ELEC_evolve::td_scf_thr*41.34);	// 2.74 is equal to E0;
					}
					else if(istep < stepcut2)
					{
						this->vext[ir] = this->vextold[ir]*2.74*cos(w_h*istep*ELEC_evolve::td_scf_thr*41.34);
					}
					else if(istep < stepcut3)
					{
						this->vext[ir] = this->vextold[ir]*2.74*(stepcut3-istep)/stepcut1*cos(w_h*istep*ELEC_evolve::td_scf_thr*41.34);
					}
				}

                //HHG of H2
				// Type 3 will be modified into more normolized type soon.
				if (ELEC_evolve::td_vexttype == 3)
				{
					const double w_h2 = 0.0428;  //a.u.
					const double w_h3 = 0.00107;  //a.u.
					const double timenow = (istep)*ELEC_evolve::td_scf_thr*41.34;
					// The parameters should be written in INPUT!

					//this->vext[ir] = this->vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow);
					//this->vext[ir] = this->vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow)*0.01944;
					this->vext[ir] = this->vextold[ir]*2.74*cos(w_h2*timenow)*sin(w_h3*timenow)*sin(w_h3*timenow);
				}

					this->vr_eff(is,ir) = this->vltot[ir] + this->vr(is, ir) + this->vext[ir];

                //std::cout << "x: " << k <<"	" << "y: " << j <<"	"<< "z: "<< i <<"	"<< "ir: " << ir << std::endl;
                //std::cout << "vext: " << this->vext[ir] << std::endl;
                //std::cout << "vrs: " << vrs(is,ir) <<std::endl;
            }
            std::cout << "vext exists" << std::endl;

            delete[] this->vextold;
            delete[] this->vext;
        }
    }

    ModuleBase::timer::tick("potential","set_vrs_tddft");
    return;
} //end subroutine set_vrs_tddft
