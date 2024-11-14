#include "hsolver.h"

#include "module_parameter/parameter.h"
namespace hsolver
{

double set_diagethr_ks(const std::string basis_type,
                       const std::string esolver_type,
                       const std::string calculation_in,
                       const std::string init_chg_in,
                       const std::string precision_flag_in,
                       const int istep,
                       const int iter,
                       const double drho,
                       const double pw_diag_thr_init,
                       const double diag_ethr_in,
                       const double nelec_in)
{
    double res_diag_ethr = diag_ethr_in;

    if (basis_type == "pw" && esolver_type == "ksdft")
    {
        // It is too complex now and should be modified.
        if (calculation_in == "nscf")
        {
            if (res_diag_ethr - 1e-2 > -1e-5)
            {
                res_diag_ethr = std::max(1e-13, 0.1 * std::min(1e-2, PARAM.inp.scf_thr / PARAM.inp.nelec));
            }
        }
        else if (iter == 1)
        {
            if (std::abs(res_diag_ethr - 1.0e-2) < 1.0e-6)
            {
                if (init_chg_in == "file" || init_chg_in == "wfc")
                {
                    //======================================================
                    // if you think that the starting potential is good
                    // do not spoil it with a louly first diagonalization:
                    // set a strict diag ethr in the input file
                    // ()diago_the_init
                    //======================================================
                    res_diag_ethr = 1.0e-5;
                }
                else
                {
                    //=======================================================
                    // starting atomic potential is probably far from scf
                    // don't waste iterations in the first diagonalization
                    //=======================================================
                    res_diag_ethr = 1.0e-2;
                }
            }

            if (calculation_in == "md" || calculation_in == "relax" || calculation_in == "cell-relax")
            {
                res_diag_ethr = std::max(res_diag_ethr, static_cast<double>(pw_diag_thr_init));
            }
        }
        else
        {
            if (iter == 2)
            {
                res_diag_ethr = 1.e-2;
            }
            res_diag_ethr = std::min(res_diag_ethr,
                                     static_cast<double>(0.1) * drho
                                         / std::max(static_cast<double>(1.0), static_cast<double>(nelec_in)));
        }

        // It is essential for single precision implementation to keep the diag ethr
        // value less or equal to the single-precision limit of convergence(0.5e-4).
        // modified by denghuilu at 2023-05-15
        if (precision_flag_in == "single")
        {
            res_diag_ethr = std::max(res_diag_ethr, static_cast<double>(0.5e-4));
        }
    }
    else
    {
        res_diag_ethr = 0.0;
    }

    return res_diag_ethr;
}


double set_diagethr_sdft(const std::string basis_type,
                         const std::string esolver_type,
                         const std::string calculation_in,
                         const std::string init_chg_in,
                         const int istep,
                         const int iter,
                         const double drho,
                         const double pw_diag_thr_init,
                         const double diag_ethr_in,
                         const int nband_in,
                         const double stoiter_ks_ne_in)
{
    double res_diag_ethr = diag_ethr_in;

    if (basis_type == "pw" && esolver_type == "sdft")
    {
        if (calculation_in == "nscf")
        {
            res_diag_ethr = std::max(std::min(1e-5, 0.1 * PARAM.inp.scf_thr / std::max(1.0, PARAM.inp.nelec)), 1e-12);
        }
        else if (iter == 1)
        {
            if (istep == 0)
            {
                if (init_chg_in == "file")
                {
                    res_diag_ethr = 1.0e-5;
                }
                res_diag_ethr = std::max(res_diag_ethr, pw_diag_thr_init);
            }
            else
            {
                res_diag_ethr = std::max(res_diag_ethr, 1.0e-5);
            }
        }
        else
        {
            if (nband_in > 0 && stoiter_ks_ne_in > 1e-6) //PARAM.inp.nbands > 0 && this->stoiter.KS_ne > 1e-6
            {
                res_diag_ethr = std::min(res_diag_ethr, 0.1 * drho / std::max(1.0, stoiter_ks_ne_in));
            }
            else
            {
                res_diag_ethr = 0.0;
            }
        }
    }
    else
    {
        res_diag_ethr = 0.0;
    }

    return res_diag_ethr;
}


double reset_diag_ethr(std::ofstream& ofs_running,
                       const std::string basis_type,
                       const std::string esolver_type,
                       const std::string precision_flag_in,
                       const double hsover_error,
                       const double drho_in,
                       const double diag_ethr_in,
                       const double nelec_in)
{

    double new_diag_ethr = 0.0;

    if (basis_type == "pw" && esolver_type == "ksdft")
    {
        ofs_running << " Notice: Threshold on eigenvalues was too large.\n";

        ModuleBase::WARNING("scf", "Threshold on eigenvalues was too large.");

        ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho_in << std::endl;
        ofs_running << " Origin diag ethr = " << diag_ethr_in << std::endl;

        new_diag_ethr = 0.1 * drho_in / nelec_in;

        // It is essential for single precision implementation to keep the diag ethr
        // value less or equal to the single-precision limit of convergence(0.5e-4).
        // modified by denghuilu at 2023-05-15
        if (precision_flag_in == "single")
        {
            new_diag_ethr = std::max(new_diag_ethr, static_cast<double>(0.5e-4));
        }
        ofs_running << " New diag ethr = " << new_diag_ethr << std::endl;
    }
    else
    {
        new_diag_ethr = 0.0;
    }

    return new_diag_ethr;
};

double cal_hsolve_error(const std::string basis_type,
                        const std::string esolver_type,
                        const double diag_ethr_in,
                        const double nelec_in)
{
    if (basis_type == "pw" && esolver_type == "ksdft")
    {
        return diag_ethr_in * static_cast<double>(std::max(1.0, nelec_in));
    }
    else
    {
        return 0.0;
    }
};

} // namespace hsolver