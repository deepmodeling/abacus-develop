#include "occupy.h"

#include "module_base/constants.h"
#include "module_base/mymath.h"
#include "module_base/parallel_reduce.h"
#include "module_elecstate/elecstate_getters.h"

Occupy::Occupy()
{
}
Occupy::~Occupy()
{
}

//===========================================================
// Four smearing methods:
// (1) do nothing
// (2) gaussian_broadening method (need 'smearing_sigma' value).
// (3) tetrahedron method
// (4) fixed_occupations
//===========================================================

bool Occupy::use_gaussian_broadening = false;
int Occupy::gaussian_type = 0;
double Occupy::gaussian_parameter = 0.01;
bool Occupy::fixed_occupations = false;

void Occupy::decision(const std::string &name, const std::string &smearing_method, const double &smearing_sigma)
{
    ModuleBase::TITLE("Occupy", "decision");
    use_gaussian_broadening = false;
    fixed_occupations = false;

    gaussian_type = 0;
    gaussian_parameter = smearing_sigma;

    if (name == "fixed")
    {
        if (gaussian_parameter != 0.0)
        {
            ModuleBase::WARNING("smearing_decision", "Fixed occupations,gauss broadening ignored");
            ModuleBase::GlobalFunc::AUTO_SET("gaussian_parameter", 0.0);
            gaussian_parameter = 0.0;
        }
    }
    else if (name == "smearing" && smearing_method == "fixed")
    {
        if (gaussian_parameter != 0.0)
        {
            ModuleBase::WARNING("smearing_decision", "Fixed occupations,gauss broadening ignored");
            ModuleBase::GlobalFunc::AUTO_SET("gaussian_parameter", 0.0);
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
        if (gaussian_parameter == 0.0)
        {
            ModuleBase::WARNING_QUIT(
                "smearing_decision",
                "Smearing requires gaussian broadening,but gaussian_parameter = 0(default value = 0.01)");
        }
        if (smearing_method == "gaussian" || smearing_method == "gauss")
        {
            gaussian_type = 0; //  0: gaussian
        }
        else if (smearing_method == "methfessel-paxton" || smearing_method == "mp")
        {
            gaussian_type = 1; // >0 Methfessel-Paxton method.
        }
        else if (smearing_method == "mp2")
        {
            gaussian_type = 2; // 2nd Methfessel-Paxton method.
        }
        else if (smearing_method == "mp3")
        {
            // acually any order Methfessel-Paxton method can be supported in Occupy::w1gauss()
            // however the parameter is string instead of int
            ModuleBase::WARNING_QUIT("occupy", "Some refactor of smearing shoule be done before supporting any order of Methfessel-Paxton method!");
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
        ModuleBase::WARNING_QUIT("occupy", "not implemented yet!");
    }
    else if (name == "from_input")
    {
        fixed_occupations = true;
    }
    else
    {
        ModuleBase::WARNING_QUIT("occupy_decision", "occupations, not implemented");
    }
    return;
}

/**
 * @brief calculates weights and fermi energy for semiconductors and insulators
 *
 * @note It is only applicable to semiconductors and insulators that have an energy gap.
 *       Here, bands are either empty or filled, and there is no electronic transfer among different k-points.
 *
 * @param nks number of k points.
 * @param wk weight of each k point (consider symmetry).
 * @param nbands number of bands.
 * @param nelec number of electrons for this spin direction.
 * @param ekb the array save the band energy.
 * @param ef output: the highest occupied Kohn-Sham level.
 * @param wg output: weight for each k, each band.
 * @param is the spin index now.
 * @param isk distinguish k point belong to which spin.
 */
void Occupy::iweights(const int nks,
                      const std::vector<double>& wk,
                      const int nbands,
                      const double& nelec,
                      const ModuleBase::matrix& ekb,
                      double& ef,
                      ModuleBase::matrix& wg,
                      const int& is, //<- is should be -1, 0, or 1. -1 means set all spins, and 0 means spin up, 1 means spin down.
                      const std::vector<int>& isk)
{
    assert(is < 2);
    double degspin = 2.0;
    if (GlobalV::NSPIN == 4)
        degspin = 1.0;
    if (is != -1)
        degspin = 1.0;

    double ib_mind = nelec / degspin;
    int ib_min = std::ceil(ib_mind);
    if (ib_min != int(ib_mind))
    {
        ModuleBase::WARNING_QUIT("iweights", "It is not a semiconductor or insulator. Please do not set 'smearing_method=fixed', and try other options.");
    }
    ef = -1e+10;

    for (int ik = 0; ik < nks; ++ik)
    {
        // when NSPIN=2, only calculate spin up or spin down with TWO_FERMI mode(nupdown != 0)
        if (GlobalV::NSPIN == 2 && isk[ik] != is && is!=-1)
        {
            continue;
        }

        for (int ib = 0; ib < nbands; ++ib)
        {
            if (ib < ib_min)
            {
                wg(ik, ib) = wk[ik];
                ef = std::max(ef, ekb(ik, ib));
            }
            else
            {
                wg(ik, ib) = 0.0;
            }
        }
    }
#ifdef __MPI
    Parallel_Reduce::gather_max_double_all(ef);
#endif

    return;
}

/**
 * @brief calculates weights with the gaussian spreading technique
 *
 * @param nks number of k points.
 * @param wk weight of each k point(symmetry considered).
 * @param nband number of bands.
 * @param nelec number of electrons.
 * @param smearing_sigma parameter input by user.
 * @param ngauss which type of smearing.
 * @param ekb band energy.
 * @param ef output: fermi level
 * @param demet output: energy correction for metal
 * @param wg output: weight of each band at each k point
 * @param is spin
 * @param isk array to point out each k belong to which spin
 */
void Occupy::gweights(const int nks,
                      const std::vector<double>& wk,
                      const int nband,
                      const double& nelec,
                      const double& smearing_sigma,
                      const int ngauss,
                      const ModuleBase::matrix& ekb,
                      double& ef,
                      double& demet,
                      ModuleBase::matrix& wg,
                      const int& is,
                      const std::vector<int>& isk)
{
    // ModuleBase::TITLE("Occupy","gweights");
    //===============================
    //  Calculate the Fermi energy ef
    //===============================
    //  call efermig
    Occupy::efermig(ekb, nband, nks, nelec, wk, smearing_sigma, ngauss, ef, is, isk);
    demet = 0.0;

    for (int ik = 0; ik < nks; ik++)
    {
        // mohan add 2011-04-03
        if (is != -1 && is != isk[ik])
            continue;

        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            //================================
            // Calculate the gaussian weights
            //================================
            // call wgauss
            wg(ik, ib) = wk[ik] * Occupy::wgauss((ef - ekb(ik, ib)) / smearing_sigma, ngauss);

            //====================================================================
            // The correct form of the band energy is  \int e n(e) de   for e<ef
            // demet is the correction to add to the sum of eigenvalues
            //====================================================================
            // Mohan fix bug 2010-1-9
            // call w1gauss
            demet += wk[ik] * smearing_sigma * Occupy::w1gauss((ef - ekb(ik, ib)) / smearing_sigma, ngauss);
        }
    }

    return;
} // end subroutine gweights

/**
 * @brief Finds the Fermi energy
 *
 * @param ekb band energy.
 * @param nband number of bands.
 * @param nks nks number of k points.
 * @param nelec number of electrons.
 * @param wk wk weight of each k point(symmetry considered).
 * @param smearing_sigma parameter input by user.
 * @param ngauss which type of smearing.
 * @param ef output: fermi level
 * @param is spin
 * @param isk array to point out each k belong to which spin
 */
void Occupy::efermig(const ModuleBase::matrix& ekb,
                     const int nband,
                     const int nks,
                     const double& nelec,
                     const std::vector<double>& wk,
                     const double& smearing_sigma,
                     const int ngauss,
                     double& ef,
                     const int& is,
                     const std::vector<int>& isk)
{
    // ModuleBase::TITLE("Occupy","efermig");
    //==================================================================
    //  Finds the Fermi energy - Gaussian Broadening (Methfessel-Paxton)
    //==================================================================
    // find bounds for the Fermi energy. Very safe choice!
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
    double elw = ekb(0, 0);
    double eup = ekb(0, nband - 1);

    for (int ik = 1; ik < nks; ik++) // do ik = 2, nks
    {
        elw = std::min(elw, ekb(ik, 0));
        eup = std::max(eup, ekb(ik, nband - 1));
    }

    eup += 2 * smearing_sigma;
    elw -= 2 * smearing_sigma;

#ifdef __MPI
    // find min and max across pools
    Parallel_Reduce::gather_max_double_all(eup);
    Parallel_Reduce::gather_min_double_all(elw);

#endif
    //=================
    // Bisection method
    //=================
    // call sumkg
    int changetime = 0;
sumkg:

    const double sumkup = Occupy::sumkg(ekb, nband, nks, wk, smearing_sigma, ngauss, eup, is, isk);
    const double sumklw = Occupy::sumkg(ekb, nband, nks, wk, smearing_sigma, ngauss, elw, is, isk);

    if (changetime > 1000)
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
        ModuleBase::WARNING_QUIT("Occupy::efermig", "ERROS in SMEARING");
    }
    else if ((sumkup - nelec) < -eps)
    {
        eup += 2 * smearing_sigma;
        ++changetime;
        goto sumkg;
    }
    else if ((sumklw - nelec) > eps)
    {
        elw -= 2 * smearing_sigma;
        ++changetime;
        goto sumkg;
    }

    for (int i = 0; i < maxiter; i++)
    {
        //======================
        // change ef value
        //======================
        ef = (eup + elw) / 2.0;
        const double sumkmid = sumkg(ekb, nband, nks, wk, smearing_sigma, ngauss, ef, is, isk);

        if (std::abs(sumkmid - nelec) < eps)
        {
            return;
        }
        else if ((sumkmid - nelec) < -eps)
        {
            elw = ef;
        }
        else
        {
            eup = ef;
        }
    }
    return;
}

/**
 * @brief This function computes the number of states under a given energy e
 *
 * @param ekb band energy.
 * @param nband number of bands.
 * @param nks number of k points.
 * @param wk weight of each k point(symmetry considered).
 * @param smearing_sigma parameter input by user
 * @param ngauss which type of smearing.
 * @param e a givern energy
 * @param is spin
 * @param isk array to point out each k belong to which spin
 * @return (double) the number of states
 */
double Occupy::sumkg(const ModuleBase::matrix& ekb,
                     const int nband,
                     const int nks,
                     const std::vector<double>& wk,
                     const double& smearing_sigma,
                     const int ngauss,
                     const double& e,
                     const int& is,
                     const std::vector<int>& isk)
{
    // ModuleBase::TITLE("Occupy","sumkg");
    double sum2 = 0.0;
    for (int ik = 0; ik < nks; ik++)
    {
        if (is != -1 && is != isk[ik])
            continue;

        double sum1 = 0.0;
        for (int ib = 0; ib < nband; ib++)
        {
            //===========================
            // call wgauss
            //===========================
            sum1 += Occupy::wgauss((e - ekb(ik, ib)) / smearing_sigma, ngauss);
        }
        sum2 += wk[ik] * sum1;
    }

    // GlobalV::ofs_running << "\n sum2 before reduce = " << sum2 << std::endl;

#ifdef __MPI
    Parallel_Reduce::reduce_double_allpool(sum2);
#endif

    // GlobalV::ofs_running << "\n sum2 after reduce = " << sum2 << std::endl;

    return sum2;
}

double Occupy::wgauss(const double &x, const int n)
{
    // ModuleBase::TITLE("Occupy","wgauss");
    //=====================================================================
    //  This function computes the approximate theta function for the
    //  iven order n, at the point x.
    //  --> (n>=0) : Methfessel-Paxton case. See PRB 40, 3616 (1989).
    //  --> (n=-1 ): Cold smearing
    //               (Marzari-Vanderbilt). See PRL 82, 3296 (1999)
    //        1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
    //  --> (n=-99): Fermi-Dirac case: 1.0/(1.0+exp(-x)).
    //=====================================================================
    //============================
    //  The value of this function
    //============================
    double wga = 0.0;
    const double maxarg = 200.0;

    //===========================
    // Fermi-Dirac(fd) smearing
    //===========================
    if (n == -99)
    {
        if (x < -maxarg)
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
    if (n == -1)
    {
        const double xp = x - 1.00 / ModuleBase::SQRT2;
        const double arg = std::min(maxarg, xp * xp);
        wga = 0.50 * erf(xp) + 1.00 / sqrt(ModuleBase::TWO_PI) * std::exp(-arg) + 0.50;
        return wga;
    }

    //====================
    // Methfessel-Paxton       //pengfei 2014-10-13
    //====================
    wga = 0.5 * (1 - erf(-x));
    // wga = gauss_freq(x * ModuleBase::SQRT2);
    //	std::cout<<"\n x="<<x<<" wga="<<wga;
    if (n == 0)
    {
        return wga;
    }

    // double hd = 0.0;

    int ni = 0;
    const double arg = std::min(maxarg, x * x);
    double hp = std::exp(-arg);
    double h0 = 1.00;
    double h1 = -2.00 * x;
    double a = 1.0 / sqrt(ModuleBase::PI);
    for (int i = 0; i < n; i++)
    {
        a = -a / (static_cast<double>(i + 1) * 4.00);
        wga = wga + a * h1 * hp;
        ++ni;
        // std::cout << " wga = " <<wga<<std::endl;
        h0 = 2.00 * (-x) * h1 - 2.00 * static_cast<double>(ni) * h0;
        ++ni;

        h1 = 2.00 * (-x) * h0 - 2.00 * static_cast<double>(ni) * h1;
    }
    return wga;
} // end function wgauss

double Occupy::w1gauss(const double &x, const int n)
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
    if (n == -99)
    {
        if (std::abs(x) <= 36.0)
        {
            const double f = 1.00 / (1.00 + exp(-x));
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

    if (n == -1)
    {
        const double xp = x - 1.00 / ModuleBase::SQRT2;
        const double arg = std::min(200.0, xp * xp);
        w1 = 1.00 / sqrt(ModuleBase::TWO_PI) * xp * std::exp(-arg);
        return w1;
    }

    //====================
    // Methfessel-Paxton
    //====================
    const double arg = std::min(200.0, x * x);
    w1 = -0.50 * std::exp(-arg) / sqrt(ModuleBase::PI);

    // std::cout << "\n x=" << x << " w1=" << w1;
    if (n == 0) // specific case : gaussian smearing.
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
    double hp = std::exp(-arg);
    double h0 = 1.00;
    double h1 = 2.00 * x;
    double a = 1.0 / sqrt(ModuleBase::PI);
    for (int i = 0; i < n; i++)
    {
        a = -a / (static_cast<double>(i + 1) * 4.00);
        ni++;
        const double h0m1 = h0;
        h0 = 2.00 * x * h1 - 2.00 * static_cast<double>(ni) * h0;
        ni++;
        h1 = 2.00 * x * h0 - 2.00 * static_cast<double>(ni) * h1;
        w1 = w1 - a * (0.50 * h0 + static_cast<double>(ni) * h0m1) * hp;
    }

    return w1;
} // end function w1gauss
