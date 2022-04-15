//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-08
//=======================

#include "read_txt_stru.h"
#include "read_txt_tools.h"

#include <algorithm>
#include <stdexcept>

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
    read_stru(const std::string &file_name)
    {
        const std::set<std::string> labels_lat = {"Lattice_Constant", "Lattice_Vectors", "Atomic_Positions"};
        std::set<std::string> labels = {"Element", "Atoms"};
        labels.insert(labels_lat.begin(), labels_lat.end());
        
        std::vector<std::vector<std::string>> stru = Read_Txt_Tools::read_file_to_vector(file_name, {"#","\\"});

        std::map<std::string, std::vector<std::string>> lat_info;
        std::map<std::string, std::map<std::string, std::vector<std::string>>> elements_info;
        std::vector<std::vector<std::string>> atoms_info;
        
        while(!stru.empty())
        {
            std::vector<std::vector<std::string>> stru_one = Read_Txt_Tools::cut_paragraph(stru, labels);
            if(stru_one[0][0]=="Element")
            {
                std::map<std::string, std::vector<std::string>> &element_info = elements_info[stru_one[0][1]];
                for(auto ptr=stru_one.begin()+1; ptr<stru_one.end(); ++ptr)
                {
                    element_info[(*ptr)[0]] = std::vector<std::string>(ptr->begin()+1, ptr->end());
                }
            }
            else if(stru_one[0][0]=="Atoms")
            {
                // [
                //  	[O, 0, 0, 0],
                //  	[H, 0, 0, 3],
                //		[Force, 1, 1, 1]
                // ]
                if(stru_one[0].size()==1)
                {
                    atoms_info.insert(atoms_info.end(), stru_one.begin()+1, stru_one.end());
                }
                else
                {
                    std::vector<std::string> atom0 = std::vector<std::string>(stru_one[0].begin()+1, stru_one[0].end());
                    atoms_info.push_back(std::move(atom0));
                    atoms_info.insert(atoms_info.end(), stru_one.begin()+1, stru_one.end());
                }
            }
            else
            {
                const auto ptr = std::find(labels_lat.begin(), labels_lat.end(), stru_one[0][0]);
                if(ptr!=labels_lat.end())
                {
                    const std::vector<std::string> data = Read_Txt_Tools::chain(stru_one);
                    lat_info[data[0]] = std::vector<std::string>(data.begin()+1, data.end());
                }
                else
                    throw std::invalid_argument(stru_one[0][0]);
            }
        }

        std::set<std::string> elements_label;
        for(const auto & i : elements_info)
            elements_label.insert(i.first);
        std::vector<std::vector<std::string>> atoms_info_new = Read_Txt_Stru::organize_atom_info(std::move(atoms_info), elements_label);	

        return make_tuple(std::move(lat_info), std::move(elements_info), std::move(atoms_info_new));
    }    

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
        const std::set<std::string> &elements_label)
    {
        std::vector<std::vector<std::string>> atoms_info_new;
        while(!atoms_info_old.empty())
        {
            const std::vector<std::vector<std::string>> atom_one = Read_Txt_Tools::cut_paragraph(atoms_info_old, elements_label);
            atoms_info_new.push_back(Read_Txt_Tools::chain(atom_one));
        }
        return atoms_info_new;
    } 
}