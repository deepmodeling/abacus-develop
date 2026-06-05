#ifndef ELEMENT_BASIS_INDEX_ORB_H
#define ELEMENT_BASIS_INDEX_ORB_H

#include "../../source_base/element_basis_index.h"
#include <vector>

	class Numerical_Orbital_Lm;
	class LCAO_Orbitals;

namespace ModuleBase
{

namespace Element_Basis_Index
{
	extern Range construct_range( const LCAO_Orbitals &orb );

	extern Range construct_range( const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb );		// orb[T][L][N]
}

}

#endif