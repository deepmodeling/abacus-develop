#ifndef EXX_ABFS_JLE_H
#define EXX_ABFS_JLE_H

#include "exx_abfs.h"

#include <vector>
#include "module_ORB/ORB_atomic_lm.h"

class Exx_Abfs::Jle
{
public:
	
	vector<
		vector<
			vector<
				Numerical_Orbital_Lm>>> jle;

	void init_jle( const double kmesh_times );

	static int Lmax;
	static double Ecut_exx;
	static double tolerence;
};

#endif	// EXX_ABFS_JLE_H
