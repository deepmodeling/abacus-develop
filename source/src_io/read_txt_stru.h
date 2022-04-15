//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-08
//=======================

#ifndef READ_TXT_STRU_H
#define READ_TXT_STRU_H

#include <tuple>
#include <vector>
#include <map>
#include <set>
#include <string>

namespace Read_Txt_Stru
{
    //lat_info:
    // {
    //  	Lattice_Constant: [30.0, Bohr],
    //  	Lattice_Vectors:
    //  		[1.0, 0.0, 0.0,
    //  		 0.0, 1.0, 0.0,
    //  		 0.0, 0.0, 1.0],
    //  	Atomic_Positions: [Cartesian, Angstrom]
    // }
    //elements_info:
    // {
    //  	O: {
    //  		Pseudo_Potantial:	[O_ONCV_PBE-1.0.upf],
    //  		Numerical_Orbital:	[orb_O.dat]
    //  	},
    //  	C: {
    //  		Pseudo_Potantial:	[C_ONCV_PBE-1.0.upf]
    //  	}
    // }
    //atoms:
    // [
    //  	[O, 0, 0, 0],
    //  	[H, 0, 0, 3, Force, 1, 1, 1]
    // ]
    std::tuple<
        std::map<std::string, std::vector<std::string>>,
        std::map<std::string, std::map<std::string, std::vector<std::string>>>,
        std::vector<std::vector<std::string>>>
    read_stru(const std::string &file_name);



    // atoms_info_old:
    //  [
    //   	[H, 0, 0, 0],
    //   	[O, 0, 0, 3],
    //      [Force, 1, 1, 1]
    //  ]
    // atoms_info_new:
    //  [
    //   	[H, 0, 0, 0],
    //   	[O, 0, 0, 3, Force, 1, 1, 1]
    //  ]
    std::vector<std::vector<std::string>> organize_atom_info(
        std::vector<std::vector<std::string>> &&atoms_info_old,
        const std::set<std::string> &elements_label);    
}

#endif