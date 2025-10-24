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

};

#endif
