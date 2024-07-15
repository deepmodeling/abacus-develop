#ifndef MD_PARA_H
#define MD_PARA_H

#include <string>

/**
 * @brief input parameters used in md
 *
 */
struct MD_para
{
    int md_nstep = 10;                 ///< md nstep
    bool md_restart = false;               ///< 1: restart MD, 0: no restart MD
    std::string md_type = "nvt";       ///< fire, nve, nvt, npt, langevin, msst
    std::string md_thermostat = "nhc"; ///< specify the thermostat: nhc, anderson, berendsen,
                                       ///< rescaling, rescale_v
    double md_dt = 1.0;                ///< Time increment (hbar/E_hartree)
    double md_tfirst = -1.0;           ///< Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    double md_tlast = -1.0;            ///< Target temperature
    int md_dumpfreq = 1;               ///< The period to dump MD information
    int md_restartfreq = 5;            ///< The period to output MD restart information
    int md_seed = -1;                  ///< random seed for MD
    int md_prec_level = 0;             ///< precision level for vc-md

    int lj_rule = 2;                     ///< combination rules used to construct the parameter matrix for LJ potential
    bool lj_eshift = false;              ///< whether to use energy shift for LJ potential
    std::vector<double> lj_rcut = {};    ///< cutoff radius of LJ potential (\AA)
    std::vector<double> lj_epsilon = {}; ///< the value of epsilon for LJ potential (eV)
    std::vector<double> lj_sigma = {};   ///< the value of sigma for LJ potential (\AA)
    std::string pot_file = "graph.pb";   ///< the filename of potential files for CMD such as DP

    int msst_direction = 2;    ///< shock direction: 0, 1, 2
    double msst_vel = 0.0;     ///< shock msst_vel (\AA/fs)
    double msst_qmass = 1.0;   ///< cell mass-like parameter (mass^2/length^4)
    double msst_vis = 0.0;     ///< artificial msst_vis (mass/length/time)
    double msst_tscale = 0.01; ///< reduction in initial temperature (0~1)

    std::string md_pmode = "iso";    ///< NPT ensemble mode: iso, aniso, tri
    std::string md_pcouple = "none"; ///< whether couple different components:
                                     ///< xyz, xy, yz, xz, none
    double md_tfreq = 0.0;           ///< Oscillation frequency, used to determine qmass
                                     ///< of thermostats coupled with particles
    double md_pfirst = -1.0;         ///< Initial pressure
    double md_plast = -1.0;          ///< Final pressure
    double md_pfreq = 0.0;           ///< Oscillation frequency, used to determine qmass
                                     ///< of thermostats coupled with barostat
    int md_tchain = 1;               ///< num of thermostats coupled with particles
    int md_pchain = 1;               ///< num of thermostats coupled with barostat

    double md_damp = 1.0; ///< Langevin damping parameter (time units)

    double md_tolerance = 100.0; ///< tolerance for velocity rescaling (K)
    int md_nraise = 1;           ///< parameters used when md_type=nvt

    bool dump_force = true;  ///< output atomic forces into the file MD_dump or
                             ///< not. liuyu 2023-03-01
    bool dump_vel = true;    ///< output atomic velocities into the file MD_dump or
                             ///< not. liuyu 2023-03-01
    bool dump_virial = true; ///< output lattice virial into the file MD_dump or
                             ///< not. liuyu 2023-03-01

    // !!!
    // They are repeated in the MD_para and Parameter, which is not good.
    double force_thr = 1.0e-3; ///< force convergence threshold in FIRE method
    bool cal_stress = false;   ///< whether calculate stress
    int my_rank = 0;           ///< MPI rank of the processor
};

#endif // MD_PARA_H