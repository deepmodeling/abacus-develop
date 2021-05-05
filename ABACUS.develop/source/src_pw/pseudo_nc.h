#ifndef PSEUDO_NC_H
#define PSEUDO_NC_H

#include "tools.h"
#include "pseudopot_upf.h"
#include "pseudo_vl.h"

class pseudo_nc: public pseudo_vl
{
	public:

	pseudo_nc();
	~pseudo_nc();

	// <PP_BETA>
	int *lll;       // lll(nbeta), angular momentum of the beta function
	int kkbeta;		// kkbeta(nbeta), point where the beta are zero

	// <PP_DIJ>
	matrix dion;	// dion(nbeta,nbeta)
	matrix betar;	//(nbeta, mesh), radial beta_{mu} functions

	// other
	int nh;         // number of beta functions per atomic type

	void set_pseudo_nc(const Pseudopot_upf &upf);
	void print_pseudo_nc(ofstream &ofs);

};

#endif // PSEUDO_NC_H
