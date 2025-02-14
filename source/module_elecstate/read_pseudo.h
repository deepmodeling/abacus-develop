#ifndef READ_PSEUDO_H
#define READ_PSEUDO_H

#include "module_cell/unitcell.h"
#include "module_cell/cal_atoms_info.h"

namespace elecstate {

    void read_pseudo(std::ofstream& ofs, UnitCell& ucell);

    // read in pseudopotential from files for each type of atom
    void read_cell_pseudopots(const std::string& fn, std::ofstream& log, UnitCell& ucell);

    void print_unitcell_pseudo(const std::string& fn, UnitCell& ucell);
    
    //===========================================
    // calculate the total number of local basis
    // Target : nwfc, lmax,
    // 			atoms[].stapos_wf
    // 			PARAM.inp.nbands
    //===========================================
    void cal_nwfc(std::ofstream& log, UnitCell& ucell,Atom* atoms);

    //======================
    // Target : meshx
    // Demand : atoms[].msh
    //======================
    void cal_meshx(int& meshx,const Atom* atoms, const int ntype);

    //=========================
    // Target : natomwfc
    // Demand : atoms[].nchi
    // 			atoms[].lchi
    // 			atoms[].oc
    // 			atoms[].na
    //=========================
    void cal_natomwfc(std::ofstream& log,int& natomwfc,const int ntype,const Atom* atoms);
}

#endif