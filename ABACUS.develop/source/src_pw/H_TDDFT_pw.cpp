#include "src_pw/potential.h"
#include "global.h"

//==========================================================
// this function aims to add external time-dependent potential 
// (eg: linear potential) used in tddft
// fuxiang add in 2017-05
//==========================================================
void Potential::set_vrs_tddft(const int istep)
{
    TITLE("Potential","set_vrs_tddft");
    timer::tick("Potential","set_vrs_tddft");

    for (int is = 0;is < NSPIN;is++)
    {
        //====================================================
        // add external linear potential, fuxiang add in 2017/05
        //====================================================

        const int td_timescale = 1;  // get the time that td_vext influences;
        if (istep >= td_timescale)
        {
            for (int i = 0;i < pw.nrxx;i++)
            {
                this->vr_eff(is, i) = this->vltot[i] + this->vr(is, i);
            }
            cout << "td_vext = 0! " << endl;
        }
        else
        {
            this->td_vextold = new double[pw.nrxx];
            this->td_vext = new double[pw.nrxx];
            const int yz = pw.ncy*pw.nczp;
            int index, i, j, k;

            for(int ir=0; ir<pw.nrxx; ++ir)
            {
                index = ir;
                i     = index / yz; // get the z, z is the fastest
                index = index - yz * i;// get (x,y)
                j     = index / pw.nczp;// get y
                k     = index - pw.nczp*j + pw.nczp_start;// get x

                if(td_vext_dire == 1)
                {
                    if (k<pw.ncx*0.05) this->td_vextold[ir] = (0.019447*k/pw.ncx-0.001069585)*ucell.lat0;
                    else if (k>=pw.ncx*0.05 && k<pw.ncx*0.95) this->td_vextold[ir] = -0.0019447*k/pw.ncx*ucell.lat0;
                    else if (k>=pw.ncx*0.95) this->td_vextold[ir] = (0.019447*(1.0*k/pw.ncx-1)-0.001069585)*ucell.lat0;
                }
                else if(td_vext_dire == 2)
                {
                    if (j<pw.ncx*0.05) this->td_vextold[ir] = (0.019447*j/pw.ncx-0.001069585)*ucell.lat0;
                    else if (j>=pw.ncx*0.05 && j<pw.ncx*0.95)	this->td_vextold[ir] = -0.0019447*j/pw.ncx*ucell.lat0;
                    else if (j>=pw.ncx*0.95) this->td_vextold[ir] = (0.019447*(1.0*j/pw.ncx-1)-0.001069585)*ucell.lat0;
                }
                else if(td_vext_dire == 3)
                {
                    if (i<pw.ncx*0.05) this->td_vextold[ir] = (0.019447*i/pw.ncx-0.001069585)*ucell.lat0;
                    else if (i>=pw.ncx*0.05 && i<pw.ncx*0.95) this->td_vextold[ir] = -0.0019447*i/pw.ncx*ucell.lat0;
                    else if (i>=pw.ncx*0.95) this->td_vextold[ir] = (0.019447*(1.0*i/pw.ncx-1)-0.001069585)*ucell.lat0;
                }

                // Gauss
/*
                const double w = 22.13;    // eV
                const double sigmasquare = 6836;
                const double timecenter = 700;
                const double timenow = (istep-timecenter)*INPUT.md_dt*41.34;
                this->td_vext[ir] = this->td_vextold[ir]*cos(w/27.2116*timenow)*exp(-timenow*timenow*0.5/(sigmasquare))*0.25;  //0.1 is modified in 2018/1/12
*/

                //HHG of H atom
/*
                if(istep < 1875)
                {
                    this->td_vext[ir] = this->td_vextold[ir]*2.74*istep/1875*cos(0.0588*istep*INPUT.md_dt*41.34);	// 2.75 is equal to E0;
                }
                else if(istep < 5625)
                {
                    this->td_vext[ir] = this->td_vextold[ir]*2.74*cos(0.0588*istep*INPUT.md_dt*41.34);
                }
                else if(istep < 7500)
                {
                    this->td_vext[ir] = this->td_vextold[ir]*2.74*(7500-istep)/1875*cos(0.0588*istep*INPUT.md_dt*41.34);
                }
*/

                //HHG of H2

                //const double timenow = (istep)*INPUT.md_dt*41.34;
                //this->td_vext[ir] = this->td_vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow);
                //this->td_vext[ir] = this->td_vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow)*0.01944;
                //this->td_vext[ir] = this->td_vextold[ir]*2.74*cos(0.0428*timenow)*sin(0.00107*timenow)*sin(0.00107*timenow);

                this->vr_eff(is,ir) = this->vltot[ir] + this->vr(is, ir) + this->td_vext[ir];

                //cout << "x: " << k <<"	" << "y: " << j <<"	"<< "z: "<< i <<"	"<< "ir: " << ir << endl;
                //cout << "td_vext: " << this->td_vext[ir] << endl;
                //cout << "vrs: " << vrs(is,ir) <<endl;
            }
            cout << "td_vext is existed!" << endl;

            delete[] this->td_vextold;
            delete[] this->td_vext;
        }
    }


    timer::tick("potential","set_vrs_tddft");
    return;
} //end subroutine set_vrs_tddft
