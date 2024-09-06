#ifndef SYSTEM_PARAMETER_H
#define SYSTEM_PARAMETER_H
#include <ctime>
#include <string>

struct System_para
{
    // ---------------------------------------------------------------
    // --------------        Parameters         ----------------------
    // ---------------------------------------------------------------
    int myrank = 0;
    int nproc = 1;
    int mypool = 0;
    int npool = 1;
    int nproc_in_pool = 1;
    std::time_t start_time = 0;

    // ---------------------------------------------------------------
    // ------------ parameters not defined in INPUT file -------------
    // ------------ but decided by INPUT parameters      -------------
    // ---------------------------------------------------------------
    bool two_fermi = false; ///< true if "nupdown" is set

    bool dos_setemin = false; ///< true: "dos_emin_ev" is set
    bool dos_setemax = false; ///< true: "dos_emax_ev" is set
    int ncx = 0, ncy = 0,
        ncz = 0;                            ///< three dimension of FFT charge/grid, same as "nx,ny,nz"
    bool out_md_control = false;            ///< true if "out_level" is set
    bool rpa_setorb = false;                ///< true if "rpa" is set
    bool gamma_only_local = false;          ///< true if "gamma_only" is true and "lcao"
                                            ///< is true; for local orbitals.
    bool double_grid = false;               ///< true if "ndx,ndy,ndz" is larger than "nx,ny,nz"
    double uramping = -10.0 / 13.6;         /// U-Ramping method (Ry)
    std::vector<double> hubbard_u = {};     ///< Hubbard Coulomb interaction parameter U (Ry)
};
#endif