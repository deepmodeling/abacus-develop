//=======================
// AUTHOR : Huanjing Gong
// DATE :   2026-03-17
//=======================

#ifndef EXXLRIDETAIL_H
#define EXXLRIDETAIL_H

#include "conv_coulomb_pot_k.h"
#include "source_base/constants.h"

#include <cmath>
#include <vector>
#include <map>
#include <string>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

	class UnitCell;
	class K_Vectors;

namespace ExxLriDetail
{
using CoulombParam
    = std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>>;

inline void trim_malloc_cache()
{
#if defined(__GLIBC__)
	malloc_trim(0);
#endif
}

extern double default_spencer_rcut(const UnitCell& ucell, const K_Vectors& kv);

extern CoulombParam build_center2_cut_coulomb_param(const CoulombParam& coulomb_param,
                                             const UnitCell& ucell,
                                             const K_Vectors& kv,
                                             bool* synthesized_rcut = nullptr);
}

#endif