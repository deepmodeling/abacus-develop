#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "SpillageStep.h"
#include "Psi_Second.h" // use to calculaet oscillation

class Metropolis
{
	public:
	Metropolis();
	~Metropolis();

	void init(void);
	
	void move_various_temperature(SpillageStep &SS);
	
	void move_one_step(const int &itemp, const int &istep, SpillageStep &Level, const int &nsteps_now);
	
	bool accept(const double &v_old, const double &v_new, const double &temperature_in);

	void update_t(void);
	void reset_temperature(const int &istep);


	// (1) start temperature
	const double& get_spi_ini_temp(void) const {return spi_ini_temp;}
	const double& get_kin_ini_temp(void) const {return kin_ini_temp;}

	// (2) cooling rate for each new temperature
	const double& get_spi_rate(void) const { return spi_rate;}
	const double& get_kin_rate(void) const { return kin_rate;}

	// (3) number of temperature 
	const int& get_spi_ntemp(void) const { return spi_ntemp;}
	const int& get_kin_ntemp(void) const { return kin_ntemp;}
	
	// (4) number of steps in each temperature
	const int& get_spi_nsteps(void) const { return spi_nsteps;}
	const int& get_kin_nsteps(void) const { return kin_nsteps;}


	Psi_Second Psi2;

	private:

	// (1)
	double spi_ini_temp; // initial temperature for spillage.
	double kin_ini_temp; // initial temperature for kinetic energy.

	// (2)
	double spi_rate; // decrease rate for spillage value
	double kin_rate; // decrease rate for kinetic energy

	// (3)
	int spi_ntemp; // number of temperatures (for spillage)
	int kin_ntemp; // number of temperatures (for kinetical energy)

	// (4)
	int spi_nsteps; // number of steps per temperature for spillage
	int kin_nsteps; // number of steps per temperature for kinetical energy


	double temperature;
	double kappa; // mohan add 2009-10-11
	double delta_kappa; // mohan add 2009-10-31
	int try_kappa;
	
	double min_spillage; // mohan add 2009-10-15
	double *kin_last; // mohan add 2009-10-15
	double *kin_temp; // mohan add 2009-10-15, temperature for each orbital.
	string states; // mohan add 2009-10-15
	int nchi; // mohan add 2009-10-15
	int output_each_istep; // mohan add 2009-10-16

	double prob_new;
	double prob_old;

	static double Boltzmann_dist( const double &value, const double &t );
	int test;
	int random_c4number;

	
};

#endif
