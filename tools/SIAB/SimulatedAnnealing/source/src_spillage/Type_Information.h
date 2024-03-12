#ifndef TYPE_INFORMATION_H
#define TYPE_INFORMATION_H
#include "common.h"

// Information of a particular type
class Type_Information
{
    public:
    Type_Information();
    ~Type_Information();
    int id;// index of element periodic table.

    int na;// number of atoms of this type.

	string *fa;// mohan add 2009-07-12
	// f: full orbtial used.
	// a: average orbital used.

	string state; // mohan add 2009-08-27, new/skip
	// state: new, start Metropolis calculation in this step.
	// state: skip, skip this Metropolis calculation,
	// used only in case RESTART=TRUE; 
	
    int lmax;// lmax used.

    int *n;// number of radial function for each L.

	int *nbase;// number of radial function of each L, sum up all the steps information before!

	int nmax;

    int nw;

	int nw2;

    void init();// allocate n, nbase, fa;
    void cal_nw();
	void cal_nmax();
};

#endif
