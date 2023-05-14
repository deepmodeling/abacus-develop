#ifndef MD_PARA_H
#define MD_PARA_H

#include <string>

/**
 * @brief input parameters used in md
 *
 */
class MD_para
{
  public:
    MD_para()
    {
        md_restart = 0;
        md_type = "nvt";
        md_thermostat = "nhc";
        md_nstep = 10;
        md_dt = 1.0;
        md_tfirst = -1.0;
        md_tlast = -1.0;
        md_dumpfreq = 1;
        md_restartfreq = 5;
        md_seed = -1;
        md_prec_level = 0;

        lj_rcut = 8.5;
        lj_epsilon = 0.01032;
        lj_sigma = 3.405;
        pot_file = "graph.pb";

        msst_direction = 2;
        msst_qmass = -1.0;
        msst_vel = 0.0;
        msst_vis = 0.0;
        msst_tscale = 0.01;

        md_pmode = "iso";
        md_pcouple = "none";
        md_tfreq = 0.0;
        md_pfirst = -1.0;
        md_plast = -1.0;
        md_pfreq = 0.0;
        md_tchain = 1;
        md_pchain = 1;

        md_damp = 1.0;

        md_tolerance = 100.0;
        md_nraise = 1;

        dump_force = true;
        dump_vel = true;
        dump_virial = true;

        force_thr = 1.0e-3;
        cal_stress = false;
        my_rank = 0;
    };
    ~MD_para(){};

    int md_nstep;              ///< md nstep
    bool md_restart;           ///< 1: restart MD, 0: no restart MD
    std::string md_type;       ///< fire, nve, nvt, npt, langevin, msst
    std::string md_thermostat; ///< specify the thermostat: nhc, anderson, berendsen, rescaling, rescale_v
    double md_dt;              ///< Time increment (hbar/E_hartree)
    double md_tfirst;          ///< Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    double md_tlast;           ///< Target temperature
    int md_dumpfreq;           ///< The period to dump MD information
    int md_restartfreq;        ///< The period to output MD restart information
    int md_seed;               ///< random seed for MD
    int md_prec_level;         ///< precision level for vc-md

    double lj_rcut;       ///< cutoff radius of LJ potential (\AA)
    double lj_epsilon;    ///< the value of epsilon for LJ potential (eV)
    double lj_sigma;      ///< the value of sigma for LJ potential (\AA)
    std::string pot_file; ///< the filename of potential files for CMD such as DP

    int msst_direction; ///< shock direction: 0, 1, 2
    double msst_vel;    ///< shock msst_vel (\AA/fs)
    double msst_qmass;  ///< cell mass-like parameter (mass^2/length^4)
    double msst_vis;    ///< artificial msst_vis (mass/length/time)
    double msst_tscale; ///< reduction in initial temperature (0~1)

    std::string md_pmode;   ///< NPT ensemble mode: iso, aniso, tri
    std::string md_pcouple; ///< whether couple different components: xyz, xy, yz, xz, none
    double md_tfreq;        ///< Oscillation frequency, used to determine qmass of thermostats coupled with particles
    double md_pfirst;       ///< Initial pressure
    double md_plast;        ///< Final pressure
    double md_pfreq;        ///< Oscillation frequency, used to determine qmass of thermostats coupled with barostat
    int md_tchain;          ///< num of thermostats coupled with particles
    int md_pchain;          ///< num of thermostats coupled with barostat

    double md_damp; ///< Langevin damping parameter (time units)

    double md_tolerance; ///< tolerance for velocity rescaling (K)
    int md_nraise;       ///< parameters used when md_type=nvt

    bool dump_force;  ///< output atomic forces into the file MD_dump or not. liuyu 2023-03-01
    bool dump_vel;    ///< output atomic velocities into the file MD_dump or not. liuyu 2023-03-01
    bool dump_virial; ///< output lattice virial into the file MD_dump or not. liuyu 2023-03-01

    double force_thr; ///< force convergence threshold in FIRE method
    bool cal_stress;  ///< whether calculate stress
    int my_rank;      ///< MPI rank of the processor
};

#endif // MD_PARA_H