#ifndef ENERGY_H
#define ENERGY_H
#include "tools.h"

class energy
{
	public:
	friend class Electrons; 
	friend class Electrons_Stochastic;//mohan add 2021-01-31 
	friend class LOOP_elec;
	friend class Local_Orbital_Charge;
	friend class Threshold_Elec;
	friend class Forces; 
	friend class Charge;
	friend class Potential;
	friend class Occupy;
	friend class wavefunc;
	friend class eximport;
	friend class Ions;
	friend class Update_input;
	friend class Force_Lo;
	friend class Exx_pw;
	friend class ON_Tests;

    energy();
    ~energy();

	void perform_dos(void);

    double etot;    	   // the total energy of the solid
    double ef;             // the fermi energy
    double ef_up;		   // spin up fermi energy
    double ef_dw;		   // spin down fermi energy
	int printe;			   // print energy every 'printe' electron iteration.
	int iter;

	public:

    // The variables needed to compute the energies
    double etot_old;
	double etot_harris;	   // Harris-Weinert-Foulkers energy
    double eband;          // the band energy
    double deband;         // correction for variational energy
	double deband_harris;
    double descf;

    double etxcc;          // the nlcc exchange and correlation
	double exx;            // the exact exchange energy.
	double evdw;			// the vdw energy				// Peize Lin add 2021.03.09

    double demet;          // correction for metals

	int out_dos;			// control dos calculation
	int out_band;                    // control band calculation  pengfei 2014-10-13

	double dos_emin_ev;
	double dos_emax_ev;
	double dos_edelta_ev;
	double bcoeff;

    double ewald();

    //=========================================================
    // delta_escf = -\int \delta rho(r) V_scf
    // Calculate the difference between Hartree and XC energy
    // at first order in the charge density difference \delta
    // rho
    //=========================================================
	void calculate_etot(void);
	void print_etot(const bool converged, const int &itep, const int &iter, 
	const double &dr2, const double &duration, const double &ethr=0, const double &avg_iter=0, bool print = true);

	void print_band(const int &ik);

	void print_format(const string &name, const double &value);

	void calculate_harris(const int &flag); //mohan add 2012-06-05

    double delta_e(void);

    void delta_escf(void);

    
    void set_exx();       // Peize Lin add 2016-12-03
    
};

#endif //energy

