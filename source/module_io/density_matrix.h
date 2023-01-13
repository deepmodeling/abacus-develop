#ifndef DENSITY_MATRIX
#define DENSITY_MATRIX

#include <string>

namespace ModuleIO
{
	void read_dm(const int &is, const std::string &fn, double*** DM, double** DM_R);
}

#endif

