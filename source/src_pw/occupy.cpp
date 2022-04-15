#include "global.h"
#include "occupy.h"
#include "../module_base/mymath.h"
#include "../src_parallel/parallel_reduce.h"

Occupy::Occupy(){}
Occupy::~Occupy(){}

//===========================================================
// Four smearing methods:
// (1) do nothing
// (2) gaussian_broadening method (need 'smearing_sigma' value).
// (3) tetrahedron method
// (4) fixed_occupations
//===========================================================

bool Occupy::use_gaussian_broadening = false;
int Occupy::gaussian_type;
double Occupy::gaussian_parameter;
bool Occupy::use_tetrahedron_method = false;
bool Occupy::fixed_occupations = false;


void Occupy::calculate_weights(void)
{
	ModuleBase::TITLE("Occupy","calculate_weights");

	// for test
//	std::cout << " gaussian_broadening = " << use_gaussian_broadening << std::endl;
//	std::cout << " tetrahedron_method = " << use_tetrahedron_method << std::endl;
//	std::cout << " fixed_occupations = " << fixed_occupations << std::endl;

	if(GlobalV::KS_SOLVER=="selinv")
	{
		GlobalV::ofs_running << " Could not calculate occupation." << std::endl;
		return;
	}

    if (	!use_gaussian_broadening &&
            !use_tetrahedron_method &&
            !fixed_occupations)
    {
        if (GlobalV::TWO_EFERMI)
        {
            iweights(GlobalC::kv.nks, GlobalC::kv.wk, GlobalV::NBANDS, GlobalC::ucell.magnet.get_nelup() , GlobalC::wf.ekb, GlobalC::en.ef_up, GlobalC::wf.wg, 0, GlobalC::kv.isk);
            iweights(GlobalC::kv.nks, GlobalC::kv.wk, GlobalV::NBANDS, GlobalC::ucell.magnet.get_neldw() , GlobalC::wf.ekb, GlobalC::en.ef_dw, GlobalC::wf.wg, 1, GlobalC::kv.isk);
			//ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
        }
        else
        {
			// -1 means don't need to consider spin.
            iweights(GlobalC::kv.nks, GlobalC::kv.wk, GlobalV::NBANDS, GlobalC::CHR.nelec, GlobalC::wf.ekb, GlobalC::en.ef, GlobalC::wf.wg, -1, GlobalC::kv.isk);
        }
    }
    else if (use_tetrahedron_method)
    {
        ModuleBase::WARNING_QUIT("calculate_weights","not implemented yet,coming soon!");
//		if(my_rank == 0)
//		{
//			tweights(GlobalC::kv.nkstot, nspin, GlobalV::NBANDS, GlobalC::CHR.nelec, ntetra,tetra, GlobalC::wf.et, GlobalC::en.ef, GlobalC::wf.wg);
//		}
    }
    else if (use_gaussian_broadening)
    {
		if (GlobalV::TWO_EFERMI)
		{
			double demet_up = 0.0;
			double demet_dw = 0.0;
			Occupy::gweights(GlobalC::kv.nks, GlobalC::kv.wk, GlobalV::NBANDS, GlobalC::ucell.magnet.get_nelup(), gaussian_parameter, gaussian_type,
			                 GlobalC::wf.ekb, GlobalC::en.ef_up, demet_up, GlobalC::wf.wg, 0, GlobalC::kv.isk);
			Occupy::gweights(GlobalC::kv.nks, GlobalC::kv.wk, GlobalV::NBANDS, GlobalC::ucell.magnet.get_neldw(), gaussian_parameter, gaussian_type,
			                 GlobalC::wf.ekb, GlobalC::en.ef_dw, demet_dw, GlobalC::wf.wg, 1, GlobalC::kv.isk);
			GlobalC::en.demet = demet_up + demet_dw;

		}
		else
		{
			// -1 means is no related to spin.
        	Occupy::gweights(GlobalC::kv.nks, GlobalC::kv.wk, GlobalV::NBANDS, GlobalC::CHR.nelec, gaussian_parameter, gaussian_type,
                 GlobalC::wf.ekb, GlobalC::en.ef, GlobalC::en.demet, GlobalC::wf.wg, -1, GlobalC::kv.isk);

		}

        //qianrui fix a bug on 2021-7-21
        Parallel_Reduce::reduce_double_allpool(GlobalC::en.demet);
    }
    else if (fixed_occupations)
    {
		// fix occupations need nelup and neldw.
		// mohan add 2011-04-03
		GlobalC::en.ef = - 1.0e+20;
		for (int ik = 0;ik < GlobalC::kv.nks ;ik++)
		{
			for (int ibnd = 0;ibnd < GlobalV::NBANDS;ibnd++)
			{
				if (GlobalC::wf.wg(ik, ibnd) > 0.0)
				{
					GlobalC::en.ef = std::max( GlobalC::en.ef, GlobalC::wf.ekb[ik][ibnd]);
				}
			}
		}
    }

    if (GlobalV::TWO_EFERMI)
    {
		Parallel_Reduce::gather_max_double_all( GlobalC::en.ef_up );
		Parallel_Reduce::gather_max_double_all( GlobalC::en.ef_dw );
    }
    else
    {
		double ebotom = GlobalC::wf.ekb[0][0];
		double etop = GlobalC::wf.ekb[0][0];
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ib++)
			{
				ebotom = min ( ebotom, GlobalC::wf.ekb[ik][ib] );
				etop = max ( etop, GlobalC::wf.ekb[ik][ib] );
			}
		}

		// parallel
		Parallel_Reduce::gather_max_double_all( GlobalC::en.ef );
		Parallel_Reduce::gather_max_double_all( etop );
		Parallel_Reduce::gather_min_double_all( ebotom );

		//not parallel yet!
//		OUT(GlobalV::ofs_running,"Top    Energy (eV)", etop * ModuleBase::Ry_to_eV);
//      OUT(GlobalV::ofs_running,"Fermi  Energy (eV)", GlobalC::en.ef * ModuleBase::Ry_to_eV);
//		OUT(GlobalV::ofs_running,"Bottom Energy (eV)", ebotom * ModuleBase::Ry_to_eV);
//		OUT(GlobalV::ofs_running,"Range  Energy (eV)", etop-ebotom * ModuleBase::Ry_to_eV);
    }

	return;
}


void Occupy::decision(const std::string &name,const std::string &smearing_method,const double &smearing_sigma)
{
	ModuleBase::TITLE("Occupy","decision");
    use_gaussian_broadening = false;
    use_tetrahedron_method = false;
    fixed_occupations = false;

    gaussian_type = 0;
    gaussian_parameter = smearing_sigma;

    if (name == "fixed")
    {
        if ( gaussian_parameter!=0.0 )
        {
            ModuleBase::WARNING("smearing_decision","Fixed occupations,gauss broadening ignored");
            ModuleBase::GlobalFunc::AUTO_SET("gaussian_parameter",0.0);
            gaussian_parameter = 0.0;
        }
    }
    else if (name == "smearing" && smearing_method == "fixed")
    {
        if ( gaussian_parameter!=0.0 )
        {
            ModuleBase::WARNING("smearing_decision","Fixed occupations,gauss broadening ignored");
            ModuleBase::GlobalFunc::AUTO_SET("gaussian_parameter",0.0);
            gaussian_parameter = 0.0;
        }
    }

	// there are four types of smearing methods:
	// (1) gaussian
	// (2) methfessel-paxton
	// (3) Marzari-Vanderbilt
	// (4) Fermi-Dirac
    else if (name == "smearing")
    {
        use_gaussian_broadening = true;
        if ( gaussian_parameter == 0.0)
        {
            ModuleBase::WARNING_QUIT("smearing_decision",
                         "Smearing requires gaussian broadening,but gaussian_parameter = 0(default value = 0.1)");
        }
        if ( smearing_method == "gaussian" || smearing_method == "gauss")
        {
            gaussian_type = 0; //  0: gaussian
        }
        else if (smearing_method == "methfessel-paxton" || smearing_method ==  "mp")
        {
            gaussian_type = 1;  // >0 Methfessel-Paxton method.
        }
        else if (smearing_method ==  "mp2")
        {
            gaussian_type = 2;  // 2nd Methfessel-Paxton method.
        }

        else if (smearing_method == "marzari-vanderbilt" || smearing_method == "cold" || smearing_method == "mv")
        {
            gaussian_type = -1;
        }
        else if (smearing_method == "fermi-dirac" || smearing_method == "fd")
        {
            gaussian_type = -99;
        }
    }
    else if (name == "tetrahedra")
    {
        use_tetrahedron_method = true;
    }
    else if (name == "from_input")
    {
        fixed_occupations = true;
    }
    else
    {
        ModuleBase::WARNING_QUIT("occupy_decision","occupations, not implemented");
    }
    return;
}


//=============================================================
// calculates weights for semiconductors and insulators
// (bands are either empty or filled)
//=============================================================
void Occupy::iweights
(
    const int nks, //number of k points.
    const std::vector<double> &wk,//weight of each k point (consider symmetry).
    const int nband,//number of bands.
    const double &nelec, //number of electrons for this spin direction.
    double **ekb,//the array save the band energy.
    double &ef,//output: the highest occupied Kohn-Sham level.
    ModuleBase::matrix &wg,//output: weight for each k, each band.
    const int &is,//the spin index now.
    const std::vector<int> &isk//distinguish k point belong to which spin.
)
{
	assert(is<2); //not include non-collinear yet!
	double degspin;
	bool conv = false;          // pengfei 2017-04-04

	degspin = (GlobalV::NSPIN == 2)? 1.0 : 2.0;
	if(GlobalV::NSPIN == 4)degspin = 1.0;//added by zhengdy-soc

	assert( degspin > 0.0);
	double ib_min = nelec/degspin;

	int ib_min1 = ( ib_min - int(ib_min) == 0) ? int(ib_min) : int(ib_min) + 1;

	for(int ik=0; ik<nks && !conv; ik++)
	{
		for(int ib=0; ib<nband && !conv; ib++)
		{
			//std::cout << " ekb=" << ekb[ik][ib] << std::endl;
			int count =0;
			ef = ekb[ik][ib];
			for(int ik1=0; ik1<nks; ik1++)
			{
				for(int ib1=0; ib1<nband; ib1++)
				{
					if( ekb[ik1][ib1] < ef || ekb[ik1][ib1] == ef )
					{
						count++;
					}
				}
			}

			if( (GlobalV::NSPIN == 2 && count == ib_min1 * nks/2)
				|| (GlobalV::NSPIN == 1 && count == ib_min1 * nks)
				|| ((GlobalV::NSPIN == 4) && count == ib_min1 * nks))
			{
				conv = true;
			}
		}
	}


	for(int ik=0; ik<nks; ik++)
	{
		for(int ib=0; ib<nband; ib++)
		{
			if (ekb[ik][ib] < ef)
			{
				wg(ik,ib) =  wk[ik];
			}
			else if(ekb[ik][ib] == ef)
			{
				wg(ik,ib) =  wk[ik] * ( ib_min + 1.0 - double(ib_min1));
			}
			else
			{
				wg(ik,ib) = 0.0;
			}
		}
	}

    if(conv == false && GlobalC::en.iter == 2)
    {
       ModuleBase::WARNING_QUIT("Occupied","not converged, change 'smearing' method.");
    }

    return;
}// end subroutine iweights


//==========================================================
// calculates weights with the gaussian spreading technique
//==========================================================
void Occupy::gweights(
    const int nks,//number of k points.
    const std::vector<double> &wk,//weight of each k point(symmetry considered).
    const int nband,//number of bands.
    const double &nelec,//number of electrons.
    const double &smearing_sigma,//parameter input by user.
    const int ngauss,//which type of smearing.
    double **ekb,//save the band energy.
    double &ef,//ouput: fermi level
    double &demet,//output: energy correction for metal
    ModuleBase::matrix &wg,//output: weight of each band at each k point.
	const int &is, // spin
	const std::vector<int> &isk) // array to point out each k belong to which spin
{
	//ModuleBase::TITLE("Occupy","gweights");
    //===============================
    // Calculate the Fermi energy ef
    //===============================
    // call efermig
    Occupy::efermig(ekb, GlobalV::NBANDS, nks, nelec, wk, smearing_sigma, ngauss, ef, is, isk);
    demet = 0.0;

    for (int ik = 0;ik < nks;ik++)
    {
		// mohan add 2011-04-03
		if(is!=-1 && is!=isk[ik]) continue;

		for (int ib = 0;ib < GlobalV::NBANDS;ib++)
		{
			//================================
			// Calculate the gaussian weights
			//================================
			// call wgauss
			wg(ik, ib) = wk [ik] * Occupy::wgauss( (ef - ekb[ik][ib] )/ smearing_sigma, ngauss);

			//====================================================================
			// The correct form of the band energy is  \int e n(e) de   for e<ef
			// demet is the correction to add to the sum of eigenvalues
			//====================================================================
			// Mohan fix bug 2010-1-9
			// call w1gauss
			demet += wk [ik] * smearing_sigma * Occupy::w1gauss((ef - ekb[ik][ib]) / smearing_sigma, ngauss);
		}
	}

    return;
} // end subroutine gweights


void Occupy::efermig
(
    double **ekb,
    const int nband,
    const int nks,
    const double &nelec,
    const std::vector<double> &wk,
    const double &smearing_sigma,
    const int ngauss,
    double &ef,
	const int &is,
	const std::vector<int> &isk
)
{
	//ModuleBase::TITLE("Occupy","efermig");
    //==================================================================
    // Finds the Fermi energy - Gaussian Broadening (Methfessel-Paxton)
    //==================================================================
    //find bounds for the Fermi energy. Very safe choice!
    const int maxiter = 300;
    const double eps = 1.0e-10;

	/*
	for(int ik=0; ik<nks; ik++)
	{
		for(int i=0; i<nband; i++)
		{
			std::cout << " ib=" << i << " ekb=" << ekb[ik][i] << std::endl;
		}
	}
	*/

	// the first 0 stands for the first k point.
    double elw = ekb[0][0];
    double eup = ekb[0][nband-1];

    for (int ik = 1;ik < nks;ik++) // do ik = 2, nks
    {
        elw = std::min(elw, ekb[ik][0]);
        eup = std::max(eup, ekb[ik][nband-1]);
    }

    eup += 2 * smearing_sigma;
    elw -= 2 * smearing_sigma; //qianrui change 2 to 5 for high temerature

#ifdef __MPI
    // find min and max across pools
	Parallel_Reduce::gather_max_double_all( eup );
	Parallel_Reduce::gather_min_double_all( elw );

#endif
    //=================
    //Bisection method
    //=================
    // call sumkg
    const double sumkup = Occupy::sumkg(ekb, nband, nks, wk, smearing_sigma, ngauss, eup, is, isk);
    const double sumklw = Occupy::sumkg(ekb, nband, nks, wk, smearing_sigma, ngauss, elw, is, isk);

    if ((sumkup - nelec) < -eps || (sumklw - nelec) > eps)
    {
		std::cout << " SOMETHING WRONG: " << std::endl;
		std::cout << " is = " << is << std::endl;
		std::cout << " eup = " << eup << std::endl;
		std::cout << " elw = " << elw << std::endl;
		std::cout << " nband = " << nband << std::endl;
		std::cout << " nelec = " << nelec << std::endl;
		std::cout << " sumkup = " << sumkup << std::endl;
		std::cout << " sumklw = " << sumklw << std::endl;
		std::cout << " sumkup - nelec = " << sumkup - nelec << std::endl;
		std::cout << " sumklw - nelec = " << sumklw - nelec << std::endl;
		ModuleBase::WARNING_QUIT("Occupy::efermig","ERROS in SMEARING");
    }

    for (int i = 0;i < maxiter;i++)
    {
        //======================
        // change ef value
        //======================
        ef = (eup + elw) / 2.0;
        const double sumkmid = sumkg(ekb, nband, nks, wk, smearing_sigma, ngauss, ef, is, isk);

        if (abs(sumkmid - nelec) < eps)
        {
            return;
        }
        else if ((sumkmid -nelec) < -eps )
        {
            elw = ef;
        }
        else
        {
            eup = ef;
        }
    }
    return;
} // end function efermig


//===================================================================
// This function computes the number of states under a given energy e
//===================================================================
double Occupy::sumkg(
	double** ekb,
	const int nband,
	const int nks,
	const std::vector<double> &wk,
	const double &smearing_sigma,
	const int ngauss,
	const double &e,
	const int &is,
	const std::vector<int> &isk
)
{
	//ModuleBase::TITLE("Occupy","sumkg");
    double sum2 = 0.0;
    for (int ik = 0;ik < nks; ik++)
	{
		if(is!=-1 && is!=isk[ik]) continue;

		double sum1 = 0.0;
		for (int ib = 0;ib < GlobalV::NBANDS; ib++)
		{
			//===========================
			// call wgauss
			//===========================
			sum1 += Occupy::wgauss(  (e - ekb[ik][ib]) / smearing_sigma, ngauss);
		}
		sum2 += wk [ik] * sum1;
	}

	//GlobalV::ofs_running << "\n sum2 before reduce = " << sum2 << std::endl;

#ifdef __MPI
	Parallel_Reduce::reduce_double_allpool(sum2);
#endif

	//GlobalV::ofs_running << "\n sum2 after reduce = " << sum2 << std::endl;

    return sum2;
} // end function sumkg


double Occupy::wgauss(const double &x,const int n)
{
	//ModuleBase::TITLE("Occupy","wgauss");
    //=====================================================================
    // This function computes the approximate theta function for the
    // iven order n, at the point x.
    // --> (n>=0) : Methfessel-Paxton case. See PRB 40, 3616 (1989).
    // --> (n=-1 ): Cold smearing
    //              (Marzari-Vanderbilt). See PRL 82, 3296 (1999)
    //       1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
    // --> (n=-99): Fermi-Dirac case: 1.0/(1.0+exp(-x)).
    //=====================================================================
    //============================
    // The value of this function
    //============================
    double wga = 0.0;
    const double maxarg = 200.0;

	//===========================
    // Fermi-Dirac(fd) smearing
    //===========================
    if (n ==  - 99)
    {
        if (x < - maxarg)
        {
            wga = 0.0;
        }
        else if (x > maxarg)
        {
            wga = 1.0;
        }
        else
        {
            wga = 1.00 / (1.0 + std::exp(-x));
        }
        return wga;
    }

    //===================
    // Cold smearing(mv)
    //===================
    if (n == - 1)
    {
        const double xp = x - 1.00 / ModuleBase::SQRT2;
        const double arg = std::min(maxarg, xp * xp);
        wga = 0.50 * erf(xp) + 1.00 / sqrt( ModuleBase::TWO_PI ) * std::exp(- arg) + 0.50;
        return wga;
    }

	//====================
    // Methfessel-Paxton       //pengfei 2014-10-13
    //====================
    wga = 0.5 * (1 - erf(-x));
    //wga = gauss_freq(x * ModuleBase::SQRT2);
    //	std::cout<<"\n x="<<x<<" wga="<<wga;
    if (n == 0)
    {
        return wga;
    }

    //double hd = 0.0;

    int ni = 0;
    const double arg = std::min(maxarg, x * x);
    double hp = std::exp(- arg);
    double h0 = 1.00;
    double h1 = -2.00 * x;
    double a = 1.0 / sqrt(ModuleBase::PI);
    for (int i = 0;i < n;i++)
    {
        a = - a / ( static_cast<double>(i+1) * 4.00);
        wga = wga + a * h1 * hp;
        ++ni;
        //std::cout << " wga = " <<wga<<std::endl;
        h0 = 2.00 * (-x) * h1 - 2.00 * static_cast<double>(ni) * h0 ;
        ++ni;

        h1 = 2.00 * (-x) * h0 - 2.00 * static_cast<double>(ni) * h1 ;
    }
    return wga;
} // end function wgauss


double Occupy::w1gauss(const double &x,const int n)
{
    //========================================================================
    //    w1gauss(x,n) = \int_{-\infty}^x   y delta(y) dy
    //    where delta(x) is the current approximation for the delta function,
    //    as obtained from w0gauss(x,n)
    //
    // --> (n>=0) : Methfessel-Paxton case
    //
    // --> (n=-1): Cold smearing (Marzari-Vanderbilt)
    //     w1gauss = 1/sqrt(2*pi)*(x-1/sqrt(2))*exp(-(x-1/sqrt(2))**2)
    //
    // --> (n=-99): Fermi-Dirac case. In this case w1gauss corresponds
    //     to the negative of the electronic entropy.
    //========================================================================
    double w1 = 0.0;

    //=======================
    // Fermi-Dirac smearing
    //=======================
    if (n ==  - 99)
    {
        if (abs(x) <= 36.0)
        {
            const double f = 1.00 / (1.00 + exp(- x));
            const double onemf = 1.00 - f;
            w1 = f * log(f) + onemf * log(onemf);
            //==================================================
            // in order to avoid problems for large values of x
            //==================================================
        }
        else
        {
            //=============================================
            // neglect w1gauss when abs(w1gauss) < 1.0d-14
            //=============================================
            w1 = 0.00;
        }
        return w1;
    }
    //===============
    // Cold smearing
    //===============

    if (n == - 1)
    {
        const double xp = x - 1.00 / ModuleBase::SQRT2;
        const double arg = std::min(200.0, xp*xp);
        w1 = 1.00 / sqrt(ModuleBase::TWO_PI) * xp * std::exp(- arg);
        return w1;
    }


	//====================
    // Methfessel-Paxton
    //====================
    const double arg = std::min(200.0, x * x);
    w1 = - 0.50 * std::exp(- arg) / sqrt(ModuleBase::PI);

	//std::cout << "\n x=" << x << " w1=" << w1;
    if (n == 0)//specific case : gaussian smearing.
    {
        return w1;
    }

/*    double hd = 0.0;
    double hp = exp(- arg);
    int ni = 0;
    double a = 1.0 / sqrt(ModuleBase::PI);

    for (int i = 0;i < n;i++)
    {
        hd = 2.00 * x * hp - 2.00 * static_cast<double>(ni) * hd;//dble(ni)
        ni ++;
        const double hpm1 = hp;
        hp = 2.00 * x * hd - 2.00 * static_cast<double>(ni) * hp;//dble(ni)
        ni ++;
		// mohan fixed bug 2010-1-10
		// description of bug: i should not be 0.
        a = - a / ( static_cast<double>(i+1) * 4.00);//dble(i)
        std::cout << " w1 == "<<w1<<std::endl;
        w1 -= a * (0.50 * hp + static_cast<double>(ni) * hpm1);//dble(ni)
        std::cout << " w1 == "<<w1<<std::endl;
    }*/

    int ni = 0;
    double hp = std::exp(- arg);
    double h0 = 1.00;
    double h1 = 2.00 * x;
    double a = 1.0 / sqrt(ModuleBase::PI);
    for (int i = 0;i < n;i++)
    {
        a = - a / ( static_cast<double>(i+1) * 4.00);
        ni++;
        const double h0m1 = h0;
        h0 =  2.00 * x * h1 - 2.00 * static_cast<double>(ni) * h0 ;
        ni++;
        h1 = 2.00 * x * h0 - 2.00 * static_cast<double>(ni) * h1 ;
        w1 = w1 - a * ( 0.50 * h0 + static_cast<double>(ni) * h0m1 ) * hp ;
    }

    return w1;
} // end function w1gauss


void Occupy::tweights(const int nks,const int nspin,const int nband,const double &nelec,
                      const int ntetra,const ModuleBase::matrix &tetra, double **ekb, double &ef, ModuleBase::matrix &wg)
{
    //===================================================================
    // calculates weights with the tetrahedron method (Bloechl version)
    // integer :: nks, nspin, GlobalV::NBANDS, ntetra, tetra (4, ntetra)
    //===================================================================

    double e1, e2, e3, e4, c1, c2, c3, c4, dosef;
    int ik, ibnd, nt, nk, ns, i, kp1, kp2, kp3, kp4;

	double etetra[4];
    int itetra[4];

    // Calculate the Fermi energy ef
    efermit(ekb, GlobalV::NBANDS, nks, nelec, nspin, ntetra, tetra, ef);

    for (ik = 0;ik < nks;ik++)
    {
        for (ibnd = 0;ibnd < nband;ibnd++)
        {
            wg(ik, ibnd) = 0.0;
        } // enddo
    } // enddo

    for (ns = 0;ns < nspin;ns++)
    {
        //==================================================================
        // nk is used to select k-points with up (ns=1) or down (ns=2) spin
        //==================================================================
        if (ns == 1)
        {
            nk = 0;
        }
        else
        {
            nk = nks / 2;
        }

        for (nt = 0;nt < ntetra;nt++)
        {
            for (ibnd = 0;ibnd < GlobalV::NBANDS;ibnd++)
            {
                //======================================================
                // etetra are the energies at the vertexes of the nt-th
                // tetrahedron
                //======================================================
                for (i = 0;i < 4;i++)
                {
                    etetra [i] = ekb[static_cast<int>( tetra(nt,i) ) + nk][ibnd];
                }

                itetra[0] = 0;

                ModuleBase::hpsort(4, etetra, itetra);

                //===============================================
                // ...sort in ascending order: e1 < e2 < e3 < e4
                //===============================================
                e1 = etetra [0];
                e2 = etetra [1];
                e3 = etetra [2];
                e4 = etetra [3];

                //==============================================================
                // kp1-kp4 are the irreducible k-points corresponding to e1-e4
                //==============================================================

                kp1 = static_cast<int>( tetra(nt,itetra[0]) )+ nk;
                kp2 = static_cast<int>( tetra(nt,itetra[1]) )+ nk;
                kp3 = static_cast<int>( tetra(nt,itetra[2]) )+ nk;
                kp4 = static_cast<int>( tetra(nt,itetra[3]) )+ nk;

                //======================
                // calculate weights wg
                //======================
                if (ef >= e4)
                {
                    wg(kp1, ibnd) = wg(kp1, ibnd) + 0.250 / ntetra;
                    wg(kp2, ibnd) = wg(kp2, ibnd) + 0.250 / ntetra;
                    wg(kp3, ibnd) = wg(kp3, ibnd) + 0.250 / ntetra;
                    wg(kp4, ibnd) = wg(kp4, ibnd) + 0.250 / ntetra;
                }
                else if (ef < e4 && ef >= e3)
                {
                    c4 = 0.250 / ntetra * pow(e4 - ef, 3) / (e4 - e1) / (e4 - e2)
                         / (e4 - e3);
                    dosef = 3.0 / ntetra * (e4 - ef) * (e4 - ef) / (e4 - e1) / (e4 - e2)
                            / (e4 - e3);
                    wg(kp1, ibnd) = wg(kp1, ibnd) + 0.250 / ntetra - c4 *
                                    (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp1][ibnd]) / 40.0;
                    wg(kp2, ibnd) = wg(kp2, ibnd) + 0.250 / ntetra - c4 *
                                    (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp2][ibnd]) / 40.0;
                    wg(kp3, ibnd) = wg(kp3, ibnd) + 0.250 / ntetra - c4 *
                                    (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp3][ibnd]) / 40.0;
                    wg(kp4, ibnd) = wg(kp4, ibnd) + 0.250 / ntetra - c4 *
                                    (4.0 - (e4 - ef) * (1.0 / (e4 - e1) + 1.0 / (e4 - e2)
                                                        + 1.0 / (e4 - e3))) + dosef * (e1 + e2 + e3 + e4 - 4.0 *
                                                                                       ekb[kp4][ibnd]) / 40.0;
                }

                else if (ef < e3 && ef >= e2)
                {
                    c1 = 0.250 / ntetra * (ef - e1) * (ef - e1) / (e4 - e1) / (e3 - e1);
                    c2 = 0.250 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef)
                         / (e4 - e1) / (e3 - e2) / (e3 - e1);
                    c3 = 0.250 / ntetra * (ef - e2) * (ef - e2) * (e4 - ef) / (e4 - e2)
                         / (e3 - e2) / (e4 - e1);
                    dosef = 1.0 / ntetra / (e3 - e1) / (e4 - e1) * (3.0 *
                            (e2 - e1) + 6.0 * (ef - e2) - 3.0 * (e3 - e1 + e4 - e2)
                            * (ef - e2) * (ef - e2) / (e3 - e2) / (e4 - e2));
                    wg(kp1, ibnd) = wg(kp1, ibnd) + c1 + (c1 + c2) * (e3 - ef)
                                    / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef *
                                    (e1 + e2 + e3 + e4 - 4.0 * ekb[kp1][ibnd]) / 40.0;
                    wg(kp2, ibnd) = wg(kp2, ibnd) + c1 + c2 + c3 + (c2 + c3)
                                    * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef *
                                    (e1 + e2 + e3 + e4 - 4.0 * ekb[kp2][ibnd]) / 40.0;
                    wg(kp3, ibnd) = wg(kp3, ibnd) + (c1 + c2) * (ef - e1)
                                    / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef *
                                    (e1 + e2 + e3 + e4 - 4.0 * ekb[kp3][ibnd]) / 40.0;
                    wg(kp4, ibnd) = wg(kp4, ibnd) + (c1 + c2 + c3) * (ef - e1)
                                    / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 +
                                            e3 + e4 - 4.0 * ekb[kp4][ibnd]) / 40.0;
                }
                else if (ef < e2 && ef >= e1)
                {
                    c4 = 0.250 / ntetra * (ef - e1) * (ef - e1) * (ef - e1) / (e2 - e1) /
                         (e3 - e1) / (e4 - e1);
                    dosef = 3.0 / ntetra * (ef - e1) * (ef - e1) / (e2 - e1) / (e3 - e1)
                            / (e4 - e1);
                    wg(kp1, ibnd) = wg(kp1, ibnd) + c4 * (4.0 - (ef - e1)
                                                          * (1.0 / (e2 - e1) + 1.0 / (e3 - e1) + 1.0 / (e4 - e1)))
                                    + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp1][ibnd]) / 40.0;
                    wg(kp2, ibnd) = wg(kp2, ibnd) + c4 * (ef - e1) / (e2 - e1)
                                    + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp2][ibnd]) / 40.0;
                    wg(kp3, ibnd) = wg(kp3, ibnd) + c4 * (ef - e1) / (e3 - e1)
                                    + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp3][ibnd]) / 40.0;
                    wg(kp4, ibnd) = wg(kp4, ibnd) + c4 * (ef - e1) / (e4 - e1)
                                    + dosef * (e1 + e2 + e3 + e4 - 4.0 * ekb[kp4][ibnd]) / 40.0;
                } // endif
            } // enddo
        } // enddo
    } // enddo

    //=====================================================================
    // add correct spin normalization : 2 for LDA, 1 for LSDA calculations
    //=====================================================================
    for (ik = 0;ik < nks;ik++)
    {
        for (ibnd = 0;ibnd < GlobalV::NBANDS;ibnd++)
        {
            wg(ik, ibnd) = wg(ik, ibnd) * 2.0 / nspin;
        }
    }
    return;
} // end subroutine tweights


double Occupy::wsweight(const ModuleBase::Vector3<double> &r, ModuleBase::Vector3<double> *rws,const int nrws)
{
    //============================================================
    // integer ir, nreq, nrws
    // real(kind=dp) r(3), rrt, ck, eps, rws(0:3,nrws), wsweight
    // parameter (eps=1.0e-6)
    //============================================================
    const double eps = 1.0e-6;

    int nreq = 1;

    for (int ir = 0;ir < nrws;ir++)
    {
        const double rrt = r * rws[ir];
        const double ck = rrt - rws[ir].x;
        //	rrt = r[1]*rws(1,ir) + r[2]*rws(2,ir) + r[3]*rws(3,ir);
        //	ck = rrt-rws(0,ir);

        if (ck > eps)
        {
            break;
        }

        if (abs(ck) < eps)
        {
            nreq++;
        }
    } // end do

    const double wswe = 1.0 / nreq;

    return wswe;
} // end function wsweight

void Occupy::efermit(double** ekb,const int nband,const int nks,const double &nelec,const int nspin,
                     const int ntetra,const ModuleBase::matrix &tetra, double &ef)
{
    //=======================================================
    // Finds the Fermi energy - tetrahedron method (Bloechl)
    // the transformation Ry to eV
    //=======================================================

    // parameter :
    const int maxiter = 300;
    const double eps = 1.0e-10;

    double efbetter;

    //===================================
    // nlw : the minimum energy band
    // elw : the lower limit of the fermi ener
    // eup : the upper limit of the fermi ener
    // external sumkt
    // find bounds for the Fermi energy.
    //===================================
    const int nlw = max(  1, static_cast<int>( (nelec / 2.0 - 5.0) )  );
    double elw = ekb[nlw][0];
    double eup = ekb[0][GlobalV::NBANDS-1];

    for (int ik = 1;ik < nks;ik++)// do ik = 2, nks
    {
        elw = min(elw, ekb[ik][nlw]);
        eup = max(eup, ekb[ik][GlobalV::NBANDS-1]);
    }
    for (int ik = 1;ik < nks;ik++)// do ik = 2, nks
    {
        elw = min(elw, ekb[ik][nlw]);
        eup = max(eup, ekb[ik][GlobalV::NBANDS-1]);
    }

    //===============================
    // Bisection method
    // the number of states with eup
    // the number of states with elw
    //===============================
    const double sumkup = sumkt(ekb, GlobalV::NBANDS, nks, nspin, ntetra, tetra, eup);
    const double sumklw = sumkt(ekb, GlobalV::NBANDS, nks, nspin, ntetra, tetra, elw);

	GlobalV::ofs_running << "\n sumkup = " << sumkup;
	GlobalV::ofs_running << "\n sumklw = " << sumklw << std::endl;

    if ((sumkup - nelec) < - eps || (sumklw - nelec) > eps)
    {
        ModuleBase::WARNING("efermit","unexpected error.");
    }

    double better = 1.0e+10;

    bool converge = false;

    double sumkmid = 0.0;
    for (int iter = 0;iter < maxiter;iter++)
    {
        // the number of states with ef
        ef = (eup + elw) / 2.0;
        sumkmid = sumkt(ekb, GlobalV::NBANDS, nks, nspin, ntetra, tetra, ef);

        if (abs(sumkmid - nelec) < better)
        {
            better = abs(sumkmid - nelec);
            efbetter = ef;
        }

        // converged
        if (abs(sumkmid - nelec) < eps)
        {
            converge = true;
            break;
        }
        else if ((sumkmid - nelec) < - eps)
        {
            elw = ef;
        }
        else
        {
            eup = ef;
        }
    }
    if (!converge)
    {
        // unconverged exit:
        // the best available ef is used . Needed in some difficult cases
        ef = efbetter;
        sumkmid = sumkt(ekb, GlobalV::NBANDS, nks, nspin, ntetra, tetra, ef);
    }

    //==============================================================
    // Check if Fermi level is above any of the highest eigenvalues
    //==============================================================
    for (int ik = 0;ik < nks;ik++)
    {
        if (ef > ekb[ik][GlobalV::NBANDS-1] + 1.e-4)
        {
            std::cout << "\n ef = " << ef;
        }
    }
    return;
} // end subroutine efermit


double Occupy::sumkt(double** ekb,const int nband,const int nks,const int nspin,const int ntetra,
                     const ModuleBase::matrix &tetra,const double &e)
{
    double etetra[4];
    double sum = 0.0;

    int nk = 0 ;
    for (int ns = 0; ns < nspin;ns++)
    {
        //==================================================================
        // nk is used to select k-points with up (ns=1) or down (ns=2) spin
        //==================================================================
        if (ns == 1)
        {
            nk = 0;
        }
        else
        {
            nk = nks / 2;
        }

        for (int nt = 0; nt < ntetra; nt++)
        {
            for (int ibnd = 0; ibnd < GlobalV::NBANDS; ibnd++)
            {
                //======================================================
                // etetra are the energies at the vertexes of the nt-th
                // tetrahedron
                //======================================================
                for (int i = 0; i < 4; i++)
                {
                    etetra [i] = ekb[  static_cast<int>( (tetra(i, nt) + nk) )][ ibnd  ];
                }

                piksort(4, etetra);
                //===========================================
                //sort in ascending order: e1 < e2 < e3 < e4
                //===========================================
                const double e1 = etetra [0];
                const double e2 = etetra [1];
                const double e3 = etetra [2];
                const double e4 = etetra [3];

                //===============================================
                // calculate sum over k of the integrated charge
                //===============================================
                if (e >= e4)
                {
                    sum += 1.0 / ntetra;
                }
                else if (e < e4 && e >= e3)
                {
                    sum += 1.0 / ntetra * (1.0 - pow((e4 - e), 3) / (e4 - e1)
                                           / (e4 - e2) / (e4 - e3));
                }
                else if (e < e3 && e >= e2)
                {
                    sum += 1.0 / ntetra / (e3 - e1) / (e4 - e1) *
                           ((e2 - e1) * (e2 - e1) + 3.0 * (e2 - e1) * (e - e2) +
                            3.0 * (e - e2) * (e - e2) - (e3 - e1 + e4 - e2) /
                            (e3 - e2) / (e4 - e2) * pow((e - e2), 3));
                }
                else if (e < e2 && e >= e1)
                {
                    sum += 1.0 / ntetra * pow((e - e1), 3) /
                           (e2 - e1) / (e3 - e1) / (e4 - e1);
                }
            }//ibnd
        }//nt
    }//ns

// add correct spin normalization : 2 for LDA, 1 for LSDA calculations
    sum *= 2.0 / nspin;
    return sum;
} // end function sumkt

void Occupy::piksort(const int n, double *a)
{
    int i;
    bool b = true;
    for (int j = 1;j < n;j++) // do j = 2, n
    {
        const double temp = a [j];
        for (i = j - 1;i >= 0;i--)  // do i = j - 1, 1, - 1
        {
            if (a [i] <= temp)
            {
                b = false;
                break;
            }
            a [i + 1] = a [i];
        }
        if (b)
        {
            i = 0;
        }
        a [i + 1] = temp;
    }
    return;
} //end subroutine piksort
