#ifndef OUT_ORBITAL_H
#define OUT_ORBITAL_H

#include "common.h"

// output orbital information
class Out_Orbital
{
	public:
	Out_Orbital();
	~Out_Orbital();

	void write(void);

	private:
	void version_information( ofstream &ofs );
	void INPUTs_information( ofstream &ofs);
	void spillage_information( ofstream &ofs, const int &ilevel);
	void metropolis_information( ofstream &ofs);
	void c4_information( ofstream &ofs);
	void mkb_information( ofstream &ofs);
	
};

#endif
