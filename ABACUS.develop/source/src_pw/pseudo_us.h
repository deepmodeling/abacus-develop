#ifndef PSEUDO_US_H
#define PSEUDO_US_H

#include "tools.h"
#include "pseudopot_upf.h"
#include "pseudo_nc.h"

class pseudo_us: public pseudo_nc
{
	public:

	pseudo_us();
	~pseudo_us();

	// <PP_QIJ>
	int nqf;

	//  <PP_RINNER>
	double *rinner;		// rinner(0:2*lmax)
	matrix qqq;			// qqq(nbeta,nbeta)
	realArray qfunc;	// qfunc(mesh,nbeta,nbeta)

	// <PP_QFCOEF>
	realArray qfcoef;	// qfcoef(nqf,0:2*lmax,nbeta,nbeta)

	// other ?
	int nqlc;			// number of angular momenta in Q

	void set_pseudo_us(const Pseudopot_upf &upf);

	void print_pseudo_us(ofstream &ofs);

};

#endif //PSEUDO_US_H
