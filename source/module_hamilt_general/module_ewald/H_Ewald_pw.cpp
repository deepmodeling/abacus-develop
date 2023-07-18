#include "H_Ewald_pw.h"
#include "module_base/mymath.h" // use heapsort
#include "dnrm2.h"
#include "module_base/parallel_reduce.h"
#include "module_base/constants.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

double H_Ewald_pw::alpha=0.0;
int H_Ewald_pw::mxr = 200;
H_Ewald_pw::H_Ewald_pw(){};
H_Ewald_pw::~H_Ewald_pw(){};

double H_Ewald_pw::compute_ewald(const UnitCell& cell,
                                 const ModulePW::PW_Basis* rho_basis,
                                 const ModuleBase::ComplexMatrix& strucFac)
{
    ModuleBase::TITLE("H_Ewald_pw","compute_ewald");
    ModuleBase::timer::tick("H_Ewald_pw","compute_ewald");

//----------------------------------------------------------
// Calculates Ewald energy with both G- and R-space terms.
// Determines optimal alpha. Should hopefully work for any structure.
//----------------------------------------------------------
    //int ng=0;
    int nr=0;
    int na=0;
    int nb=0;
    //int nt=0;
    int nrm=0;

    double ewaldg=0.0;
    double ewaldr=0.0;
    double ewalds=0.0;

    ModuleBase::Vector3<double> dtau ;
    ModuleBase::Vector3<double> *r;
    double *r2;
    double rmax=0.0;
    double rr=0.0;
    double upperbound=0.0;
    double fact=0;
    // total ionic charge in the cell
    // ewald energy computed in reciprocal space
    // ewald energy computed in real space
    // the difference tau_s - tau_s'
    // alpha term in ewald sum
    // input of the rgen routine ( not used here )
    // the square modulus of R_j-tau_s-tau_s'
    // the maximum radius to consider real space sum
    // buffer variable
    // used to optimize alpha

	if(GlobalV::test_energy)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"mxr",mxr);
    r  = new ModuleBase::Vector3<double>[mxr];
    r2 = new double[mxr];
    int* irr = new int[mxr];

    // (1) calculate total ionic charge
    double charge = 0.0;
    for (int it = 0;it < cell.ntype;it++)
    {
        charge += cell.atoms[it].na * cell.atoms[it].ncpp.zv;//mohan modify 2007-11-7
    }
    if(GlobalV::test_energy)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total ionic charge",charge);

	// (2) calculate the converged value: alpha
    H_Ewald_pw::alpha = 2.90;
    do
    {
        alpha -= 0.10;
        // choose alpha in order to have convergence in the sum over G
        // upperbound is a safe upper bound for the error in the sum over G

        if (alpha <= 0.0)
        {
            ModuleBase::WARNING_QUIT("ewald","Can't find optimal alpha.");
        }
        upperbound = 2.0 * charge * charge * sqrt(2.0 * alpha / ModuleBase::TWO_PI) *
                     erfc(sqrt(cell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    }
    while (upperbound > 1.0e-7);
    if(GlobalV::test_energy)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"alpha",alpha);
	if(GlobalV::test_energy)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Upper bound",upperbound);

    // G-space sum here.
    // Determine if this processor contains G=0 and set the constant term
    if (rho_basis->ig_gge0 >= 0)
    {
        ewaldg = - charge * charge / alpha / 4.0;
    }
    else
    {
        ewaldg = 0.0;
    }

	// in plane wave basis, only use k=0 point is not 
	// called "gamma_only", only if the wave functions
	// are stored as double type, the gamma_only = true. 
	// I don't know why "gamma_only" in plane wave 
	// makes the fact below is 2, that's a little complicated
	// to understand. I think that may because only half
	// the G vectors are used. Unfortunately implement the 
	// function hasn't in my plan list yet.
	//
	// but that's not the term "gamma_only" I want to use in LCAO,  
	fact = 1.0;

    //GlobalV::ofs_running << "\n pwb.gstart = " << pwb.gstart << std::endl;

    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if(ig == rho_basis->ig_gge0) continue;
        std::complex<double> rhon = ModuleBase::ZERO;
        for (int it=0; it<cell.ntype; it++)
        {
            rhon += static_cast<double>(cell.atoms[it].ncpp.zv) * conj(strucFac(it, ig));
        }
        ewaldg += fact * std::abs(rhon) * std::abs(rhon)
                  * exp(- rho_basis->gg[ig] * cell.tpiba2 / alpha / 4.0 ) / rho_basis->gg[ig] / cell.tpiba2;
    }

    ewaldg = ModuleBase::FOUR_PI / cell.omega * ewaldg;

//	std::cout << "\n ewaldg = " << ewaldg;

    //  Here add the other constant term
	if (rho_basis->ig_gge0 >= 0)
	{
    	for (int it = 0; it < cell.ntype;it++)
    	{
        	ewaldg = ewaldg - cell.atoms[it].na * cell.atoms[it].ncpp.zv * cell.atoms[it].ncpp.zv * sqrt(8.0 / ModuleBase::TWO_PI * alpha);
		}
    }//mohan modify 2007-11-7, 2010-07-26

    // R-space sum here (only done for the processor that contains G=0)
    ewaldr = 0.0;
    if (rho_basis->ig_gge0 >= 0)
    {	
        rmax = 4.0 / sqrt(alpha) / cell.lat0;
		if(GlobalV::test_energy)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"rmax(unit lat0)",rmax);
        // with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
        int nt1=0;
        int nt2=0;

        for (nt1 = 0; nt1 < cell.ntype; nt1++)
        {
            for (nt2 = 0; nt2 < cell.ntype; nt2++)
            {
                for (na = 0; na < cell.atoms[nt1].na; na++)
                {
                    for (nb = 0; nb < cell.atoms[nt2].na; nb++)
                    {
                        //calculate tau[na]-tau[nb]
                        dtau = cell.atoms[nt1].tau[na] - cell.atoms[nt2].tau[nb];
                        //generates nearest-neighbors shells
                        H_Ewald_pw::rgen(dtau, rmax, irr, cell.latvec, cell.G, r, r2, nrm);
                        // at-->cell.latvec, bg-->G
                        // and sum to the real space part

                        if (GlobalV::test_energy>1)
                        {
                            ModuleBase::GlobalFunc::OUT("dtau.x",dtau.x);
                            ModuleBase::GlobalFunc::OUT("dtau.y",dtau.y);
                            ModuleBase::GlobalFunc::OUT("dtau.z",dtau.z);
                            ModuleBase::GlobalFunc::OUT("nrm",nrm);
                        }
                        for (nr = 0;nr < nrm;nr++)
                        {
                            rr = sqrt(r2 [nr]) * cell.lat0;
                            ewaldr = ewaldr + cell.atoms[nt1].ncpp.zv * cell.atoms[nt2].ncpp.zv *
                                     erfc(sqrt(alpha) * rr) / rr;

                        } // enddo
                        if (GlobalV::test_energy>1) ModuleBase::GlobalFunc::OUT("ewaldr",ewaldr);
                    } // enddo
                } // enddo
            } // nt2
        }//nt1
    } // endif

    ewalds = 0.50 * ModuleBase::e2 * (ewaldg + ewaldr);

	// mohan fix bug 2010-07-26
    Parallel_Reduce::reduce_double_pool( ewalds );

    if (GlobalV::test_energy>1)
    {
        ModuleBase::GlobalFunc::OUT("ewaldg",ewaldg);
        ModuleBase::GlobalFunc::OUT("ewaldr",ewaldr);
        ModuleBase::GlobalFunc::OUT("ewalds",ewalds);
    }

    delete[] irr;
    delete[] r;
    delete[] r2;

    ModuleBase::timer::tick("H_Ewald_pw","compute_ewald");
    return ewalds;
} // end function ewald


void H_Ewald_pw::rgen(
    const ModuleBase::Vector3<double> &dtau,
    const double &rmax,
    int *irr,
    const ModuleBase::Matrix3 &latvec,
    const ModuleBase::Matrix3 &G,
    ModuleBase::Vector3<double> *r,
    double *r2,
    int &nrm)
{
    //-------------------------------------------------------------------
    // generates neighbours shells (in units of alat) with length
    // less than rmax,and returns them in order of increasing length.
    // r=i*a1+j*a2+k*a3-dtau,
    // where a1,a2,a3 are the vectors defining the lattice

    // first the dummy variables
    // int nrm, mxr;
    // output: the number of vectors in the spher
    // input: the maximum number of vectors

    // real r (3, mxr), r2 (mxr), at (3, 3), bg (3, 3), dtau (3),rmax;
    // output: coordinates of vectors R+tau_s-tau
    // output: square modulus of vectors R+tau_s-
    // input: direct lattice vectors
    // input: reciprocal lattice vectors
    // input: the std::vector tau_s-tau_s'
    // input: the radius of the sphere in real sp
    //    and here the local variables

    int nm1=0;
    int nm2=0;
    int nm3=0;
    int i=0;
    int j=0;
    int k=0;
    // index on R vectors for order
    //  maximum values for trial vectors
    //  counters on trial vectors
    // counter on polarizations
    // counter on R vectors
    // index of swapping
    // used for swapping

    ModuleBase::Vector3<double> t;
    ModuleBase::Vector3<double> t1;
    double tt=0.0;
    double bg1[3]={0,0,0};
    // buffer contains the actual r
    // buffer cotains the modulus of actual r
    // used for swapping
    // function to find the norm of a std::vector
    // external dnrm2

    nrm = 0;

    if (rmax == 0.0)
    {
        return;
    }

    bg1[0] = G.e11;
    bg1[1] = G.e12;
    bg1[2] = G.e13;

    nm1 = (int)(dnrm2(3, bg1, 1) * rmax + 2);

    bg1[0] = G.e21;
    bg1[1] = G.e22;
    bg1[2] = G.e23;

    nm2 = (int)(dnrm2(3, bg1, 1) * rmax + 2);

    bg1[0] = G.e31;
    bg1[1] = G.e32;
    bg1[2] = G.e33;

    nm3 = (int)(dnrm2(3, bg1, 1) * rmax + 2);

    if (GlobalV::test_energy>1)
    {
        ModuleBase::GlobalFunc::OUT("nm1",nm1);
        ModuleBase::GlobalFunc::OUT("nm2",nm2);
        ModuleBase::GlobalFunc::OUT("nm3",nm3);
    }

    for (i = -nm1; i <= nm1; i++) // mohan fix bug, add '='. 2009-02-27
    {
        for (j = -nm2; j <= nm2; j++)
        {
            for (k = -nm3; k <= nm3; k++)
            {
                ModuleBase::Vector3<double> t1(i,j,k);
//				out.printV3(t1);
                t = t1 * latvec; // bug ! first '*latvec', second '-dtau'.
                t = t - dtau; // bug ! t = t - dtau, not t1 = t1 -tau;

//				out.printV3(t);	// mohan fix 2bugs here, 2009-2-27
//				out.printM3("latvec",latvec);

                tt = t.x * t.x + t.y * t.y + t.z * t.z;

                if (tt <= rmax * rmax && std::abs(tt) > 1.e-10)
                {
                    if (nrm > mxr)
                    {
                        std::cerr << "\n rgen, too many r-vectors," << nrm;
                    }
                    r[nrm] = t;
                    r2[nrm] = tt;
                    nrm++;
                } // endif
            } // enddo
        } // enddo
    } // enddo

    //   reorder the vectors in order of increasing magnitude
    //   initialize the index inside sorting routine
    irr[0] = 0;
    if (nrm > 1)
    {
        ModuleBase::heapsort(nrm, r2, irr);
    }

	// mohan fix bug 2011-06-07
	for(int i=0; i<nrm; i++)
	{
		back:
		const int index = irr[i];
		// large index before 'i' should be eliminated.
		// according to change elements again and again.
		if( index < i )
		{
			double swap_x = r[index].x;	
			double swap_y = r[index].y;	
			double swap_z = r[index].z;
			
			r[index].x = r[ irr[index] ].x;	
			r[index].y = r[ irr[index] ].y;			
			r[index].z = r[ irr[index] ].z;			
				
			r[ irr[index] ].x = swap_x;
			r[ irr[index] ].y = swap_y;
			r[ irr[index] ].z = swap_z;


			const int iswap = irr[i];
			irr[i] = irr[index];
			irr[index] = iswap;

			goto back;
		}	
	}

    return;
} //end subroutine rgen
