#ifndef MD_PARAMETERS_H
#define MD_PARAMETERS_H

#include "../module_base/constants.h"

class MD_parameters
{
public:
    MD_parameters()
	{
		md_restart = 0;
		md_type = 1;
        md_nstep = 10;
		md_dt = 1;
		md_tfirst = 0;
		md_tlast = -1;
		md_dumpfreq = 1;
		md_restartfreq = 1;

		// Classic MD 
		md_ensolver = "FP";
		lj_rcut = 8.5;        
		lj_epsilon = 0.01032; 
		lj_sigma = 3.405;     

		// MSST
		msst_direction = 2;
		msst_qmass = 1;
		msst_vel = 0;
		msst_vis = 0;
		msst_tscale = 0.01;

		// NHC
		md_tfreq = 1;
		md_mnhc = 4;

		// Langevin
		md_damp = 1;
	};
    ~MD_parameters(){};

    int md_nstep;                 // md nstep
    bool md_restart;              // 1: restart MD, 0: no restart MD
	int md_type;                  // -1: FIRE, 0: NVE, 1: NVT NHC, 2: LGV, 3: NVT ADS, 4: MSST 
    double md_dt;                 // Time increment (hbar/E_hartree)
    double md_tfirst;             // Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    double md_tlast;              // Target temperature
    int md_dumpfreq;              // The period to dump MD information
	int md_restartfreq;           // The period to output MD restart information

	// Classic MD               // liuyu 2021-07-30
	std::string md_ensolver;    // choose potential: LJ, DP, FP
	double lj_rcut;             // cutoff radius of LJ potential (\AA)
	double lj_epsilon;          // the value of epsilon for LJ potential (eV)
	double lj_sigma;            // the value of sigma for LJ potential (\AA)

	// MSST
	int msst_direction;              // shock direction: 0, 1, 2
	double msst_vel;            // shock msst_vel (\AA/fs)
	double msst_qmass;               // cell mass-like parameter (mass^2/length^4)
	double msst_vis;           // artificial msst_vis (mass/length/time)
	double msst_tscale;              // reduction in initial temperature (0~1)

	// NHC
	double md_tfreq;               // Oscillation frequency, used to determine msst_qmass of NHC
	int md_mnhc;                   // num of NHC

	// Langevin
	double md_damp;                // damping parameter (time units)
};


#endif