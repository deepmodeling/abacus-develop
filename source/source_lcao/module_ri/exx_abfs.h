#ifndef EXX_ABFS_H
#define EXX_ABFS_H

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <string>

#include "../../source_basis/module_ao/ORB_atomic_lm.h"
#include "../../source_base/element_basis_index.h"
#include "../../source_base/matrix.h"
#include "../../source_base/vector3.h"

class Exx_Abfs
{
public:
	class Jle;
	class IO;
	class Construct_Orbs;
	class PCA;
	
	int rmesh_times = 5;				// Peize Lin test
	int kmesh_times = 1;				// Peize Lin test

	static int get_Lmax(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb)
	{
		int Lmax = -1;
		for( const auto &orb_T : orb )
			{ Lmax = std::max( Lmax, static_cast<int>(orb_T.size())-1 ); }
		return Lmax;
	}
};

#endif
