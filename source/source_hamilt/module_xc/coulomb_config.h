#ifndef COULOMB_CONFIG_H
#define COULOMB_CONFIG_H

#include <vector>
#include <map>
#include <string>

namespace Conv_Coulomb_Pot_K
{
    enum class Coulomb_Type{Fock, Erfc};
    enum class Ccp_Type{		//	parameter:
        Ccp,					//
        Hf,						//		"hf_Rcut"
        Erfc,					//		"hse_omega"
        Erf};					//		"hse_omega", "hf_Rcut"
    enum class Coulomb_Method{Center2, Ewald}; // Different methods for constructing the Coulomb matrix.
}

using CoulombParam = std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>>;

#endif // COULOMB_CONFIG_H
