#ifndef EXX_ABFS_H
#define EXX_ABFS_H

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <string>

#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include "../module_base/element_basis_index.h"
#include "../module_base/matrix.h"
#include "../module_base/vector3.h"

class Exx_Abfs
{
public:
	class Util;
	class Abfs_Index;
	class Jle;
	class Inverse_Matrix_Double;
	class IO;
	class Construct_Orbs;
	class PCA;

	class Matrix_Orbs11;
	class Matrix_Orbs21;
	class Matrix_Orbs22;
	
	int rmesh_times = 5;				// Peize Lin test
	int kmesh_times = 1;				// Peize Lin test

};

#endif
