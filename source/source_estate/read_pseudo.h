#ifndef READ_PSEUDO_H
#define READ_PSEUDO_H

#include "source_cell/unitcell.h"

namespace elecstate {

    void read_pseudo(std::ofstream& ofs, UnitCell& ucell,
                   const std::string& pseudo_dir,
                   const std::string& global_out_dir,
                   const bool out_element_info,
                   const std::string& dft_functional,
                   const bool lspinorb,
                   const double pseudo_rcut,
                   const double soc_lambda,
                   const int nspin,
                   const int npol,
                   const std::string& basis_type,
                   const std::string& esolver_type,
                   const std::string& init_wfc,
                   const int nbands);

    // read in pseudopotential from files for each type of atom
    void read_cell_pseudopots(const std::string& fn, std::ofstream& log, UnitCell& ucell,
                              const std::string& global_out_dir,
                              const std::string& dft_functional,
                              const bool lspinorb,
                              const double pseudo_rcut,
                              const double soc_lambda);

    void print_unitcell_pseudo(const std::string& fn, UnitCell& ucell);
    
    //===========================================
    // calculate the total number of local basis
    // Target : nwfc, lmax,
    // 			atoms[].stapos_wf
    // 			PARAM.inp.nbands
    //===========================================
    void cal_nwfc(std::ofstream& log, UnitCell& ucell,Atom* atoms, const int nspin, const int nlocal, const int npol,
               const std::string& basis_type, const std::string& esolver_type, const std::string& init_wfc, const int nbands);

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
    void cal_natomwfc(std::ofstream& log,int& natomwfc,const int ntype,const Atom* atoms,const int nspin);
}

#endif