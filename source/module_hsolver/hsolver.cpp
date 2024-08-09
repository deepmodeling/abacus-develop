#include "hsolver.h"

namespace hsolver
{

double reset_diag_ethr(std::ofstream& ofs_running,
                       const double hsover_error,
                       const double drho,
                       double diag_ethr_in,
                       std::string basis_type,
                       std::string esolver_type)
{

    if (basis_type == "pw" && esolver_type == "ksdft")
    {
        ofs_running << " Notice: Threshold on eigenvalues was too large.\n";

        ModuleBase::WARNING("scf", "Threshold on eigenvalues was too large.");

        ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho << std::endl;
        ofs_running << " Origin diag ethr = " << diag_ethr_in << std::endl;

        diag_ethr_in = 0.1 * drho / GlobalV::nelec;

        // It is essential for single precision implementation to keep the diag ethr
        // value less or equal to the single-precision limit of convergence(0.5e-4).
        // modified by denghuilu at 2023-05-15
        if (GlobalV::precision_flag == "single")
        {
            diag_ethr_in = std::max(diag_ethr_in, static_cast<double>(0.5e-4));
        }
        ofs_running << " New diag ethr = " << diag_ethr_in << std::endl;
    }
    else
    {
        diag_ethr_in = 0.0;
    }

    return diag_ethr_in;
};

double cal_hsolve_error(const double diag_ethr_in, std::string basis_type, std::string esolver_type)
{
    if (basis_type == "pw" && esolver_type == "ksdft")
    {
        return diag_ethr_in * static_cast<double>(std::max(1.0, GlobalV::nelec));
    }
    else
    {
        return 0.0;
    }
};


// double set_diagethr(double diag_ethr_in,
//                   const int istep,
//                   const int iter,
//                   const double drho,
//                   std::string basis_type,
//                   std::string esolver_type)
// {
//     if (basis_type = "pw" && esolver_type = "ksdft")
//     {
//         // It is too complex now and should be modified.
//         if (iter == 1)
//         {
//             if (std::abs(diag_ethr_in - 1.0e-2) < 1.0e-6)
//             {
//                 if (GlobalV::init_chg == "file")
//                 {
//                     //======================================================
//                     // if you think that the starting potential is good
//                     // do not spoil it with a louly first diagonalization:
//                     // set a strict diag ethr in the input file
//                     // ()diago_the_init
//                     //======================================================
//                     diag_ethr_in = 1.0e-5;
//                 }
//                 else
//                 {
//                     //=======================================================
//                     // starting atomic potential is probably far from scf
//                     // don't waste iterations in the first diagonalization
//                     //=======================================================
//                     diag_ethr_in = 1.0e-2;
//                 }
//             }

//             if (GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
//             {
//                 diag_ethr_in = std::max(diag_ethr_in, static_cast<double>(GlobalV::PW_DIAG_THR));
//             }
//         }
//         else
//         {
//             if (iter == 2)
//             {
//                 diag_ethr_in = 1.e-2;
//             }
//             diag_ethr_in = std::min(diag_ethr_in,
//                                     static_cast<double>(0.1) * drho
//                                         / std::max(static_cast<double>(1.0), static_cast<double>(GlobalV::nelec)));
//         }
//         // It is essential for single precision implementation to keep the diag ethr
//         // value less or equal to the single-precision limit of convergence(0.5e-4).
//         // modified by denghuilu at 2023-05-15
//         if (GlobalV::precision_flag == "single")
//         {
//             diag_ethr_in = std::max(diag_ethr_in, static_cast<double>(0.5e-4));
//         }
//     }
//     else if (basis_type = "pw" && esolver_type = "sdft")
//     {
//         if (iter == 1)
//         {
//             if (istep == 0)
//             {
//                 if (GlobalV::init_chg == "file")
//                 {
//                     diag_ethr_in = 1.0e-5;
//                 }
//                 diag_ethr_in = std::max(diag_ethr_in, GlobalV::PW_DIAG_THR);
//             }
//             else
//             {
//                 diag_ethr_in = std::max(diag_ethr_in, 1.0e-5);
//             }
//         }
//         else
//         {
//             if (GlobalV::NBANDS > 0 && this->stoiter.KS_ne > 1e-6)
//             {
//                 diag_ethr_in = std::min(diag_ethr_in, 0.1 * drho / std::max(1.0, this->stoiter.KS_ne));
//             }
//             else
//             {
//                 diag_ethr_in = 0.0;
//             }
//         }
//     }
//     else
//     {
//         diag_ethr_in = 0.0;
//     }
     
//     return 0.0;
// };


} // namespace hsolver