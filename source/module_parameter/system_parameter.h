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
    int nthread_per_proc = 1;
    int mypool = 0;
    int npool = 1;
    int nproc_in_pool = 1;
    std::time_t start_time = 0;

    // ---------------------------------------------------------------
    // ------------ parameters not defined in INPUT file -------------
    // ------------ but decided by INPUT parameters      -------------
    // ---------------------------------------------------------------
    int nlocal = 0; ///< total number of local basis.
    bool two_fermi = false; ///< true if "nupdown" is set
    bool use_uspp = false;   ///< true if "uspp" is set
    bool dos_setemin = false; ///< true: "dos_emin_ev" is set
    bool dos_setemax = false; ///< true: "dos_emax_ev" is set

    double dq = 0.010; // space between Q points of the reciprocal radial tab
    int nqx = 10000;   // number of points describing reciprocal radial tab
    int nqxq = 10000;  // number of points describing reciprocal radial tab for Q

    int ncx = 0, ncy = 0,
        ncz = 0;                            ///< three dimension of FFT charge/grid, same as "nx,ny,nz"
    bool out_md_control = false;            ///< true if "out_level" is set
    bool gamma_only_pw = false;             ///< true if "gamma_only" is true and "basis_type" is "pw"
                                            ///< for plane wave basis.
    bool gamma_only_local = false;          ///< true if "gamma_only" is true and "lcao"
                                            ///< is true; for local orbitals.
    std::string global_in_card = "INPUT";   ///< input file
    std::string global_in_stru = "STRU";    ///< stru file
    std::string global_out_dir = "";        ///< global output directory
    std::string global_readin_dir = "";     ///< global readin directory
    std::string global_stru_dir = "";       ///< global structure directory
    std::string global_matrix_dir = "";     ///< global matrix directory
    std::string log_file = "log";           ///< log file name

    bool deepks_setorb = false;             ///< true if "deepks" is set
    int npol = 1;                           ///< number of polarization
    bool domag = false;                     /// 1 : calculate the magnetism with x, y, z component
    bool domag_z = false;                   /// 1 : constrain the magnetism to z axis

    bool double_grid = false;               ///< true if "ndx,ndy,ndz" is larger than "nx,ny,nz"
    double uramping = -10.0 / 13.6;         /// U-Ramping method (Ry)
    std::vector<double> hubbard_u = {};     ///< Hubbard Coulomb interaction parameter U (Ry)
    int kpar_lcao = 1;                      ///< global number of pools for LCAO diagonalization only
    int nbands_l = 0;                       ///< number of bands of each band parallel calculation, same to nbands when bndpar=1
    bool ks_run = false;                    ///< true if current process runs KS calculation
    bool all_ks_run = true;                 ///< true if only all processes run KS calculation
};
#endif