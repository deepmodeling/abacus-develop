#ifndef MULTI_ZETA_H
#define MULTI_ZETA_H

#include "common.h"
#include "SpillageStep.h"
#include "Metropolis.h"
#include "Orthogonal.h"

class MultiZeta
{
	public:
	MultiZeta();
	~MultiZeta();

	SpillageStep *Level;
	Metropolis metro;
	Orthogonal Orth; //Orthogonal wave functions.

	void init(void);

	int nlevel;

	int *lmax_type; // lmax for each type(consider all levels)

	IntArray l_nchi;

	// mohan add 2010-05-02
	// convinent for output information
	int ilevel;


	void set_levels(void);
	void cal_lmax_nmax(void);

	private:
	int test;

	// save the coefficients after each level.
	// mohan add 2010-04-11
	void saveC(void);

};

#endif
