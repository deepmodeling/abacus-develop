/**
 * @file read_pp_ucell.h
 * @brief Functions for reading pseudopotential files.
 */
#ifndef READ_PP_UCELL_H
#define READ_PP_UCELL_H

#include "unitcell.h"
#include "cal_atoms_info.h"

namespace unitcell {

    /**
     * @brief Read pseudopotential files and calculate atom information.
     *
     * @param ofs output file stream [in]
     * @param ucell unit cell [in/out]
     * @param pseudo_dir pseudopotential directory [in]
     * @param global_out_dir global output directory [in]
     * @param out_element_info whether to output element information [in]
     * @param dft_functional DFT functional [in]
     * @param lspinorb spin-orbit coupling flag [in]
     * @param pseudo_rcut pseudo cut-off radius [in]
     * @param soc_lambda SOC lambda parameter [in]
     * @param nspin number of spin components [in]
     * @param npol number of polarizations [in]
     * @param basis_type basis type [in]
     * @param esolver_type solver type [in]
     * @param init_wfc initial wavefunction type [in]
     * @param nbands number of bands [in]
     * @param two_fermi two Fermi levels flag [in]
     * @param nelec_delta electron number delta [in]
     * @param smearing_method smearing method [in]
     * @param ks_solver KS solver type [in]
     * @param bndpar band parallel parameter [in]
     * @param nelec number of electrons [in]
     * @param nupdown spin polarization [in]
     * @return AtomsInfoResult containing calculated atom information
     */
    AtomsInfoResult read_pseudo(std::ofstream& ofs, UnitCell& ucell,
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
                   const int nbands,
                   const bool two_fermi,
                   const double nelec_delta,
                   const std::string& smearing_method,
                   const std::string& ks_solver,
                   const int bndpar,
                   const double nelec,
                   const double nupdown);

    /**
     * @brief Read pseudopotential from files for each type of atom.
     *
     * @param fn filename [in]
     * @param log output file stream [in]
     * @param ucell unit cell [in/out]
     * @param global_out_dir global output directory [in]
     * @param dft_functional DFT functional [in]
     * @param lspinorb spin-orbit coupling flag [in]
     * @param pseudo_rcut pseudo cut-off radius [in]
     * @param soc_lambda SOC lambda parameter [in]
     */
    void read_cell_pseudopots(const std::string& fn, std::ofstream& log, UnitCell& ucell,
                              const std::string& global_out_dir,
                              const std::string& dft_functional,
                              const bool lspinorb,
                              const double pseudo_rcut,
                              const double soc_lambda);

    /**
     * @brief Print unit cell pseudopotential information.
     *
     * @param fn filename [in]
     * @param ucell unit cell [in]
     */
    void print_unitcell_pseudo(const std::string& fn, UnitCell& ucell);
    
    /**
     * @brief Calculate the total number of local basis.
     *
     * Target: nwfc, lmax, atoms[].stapos_wf, PARAM.inp.nbands
     *
     * @param log output file stream [in]
     * @param ucell unit cell [in/out]
     * @param atoms atom pointer [in/out]
     * @param nspin number of spin components [in]
     * @param nlocal total number of local basis [in]
     * @param npol number of polarizations [in]
     * @param basis_type basis type [in]
     * @param esolver_type solver type [in]
     * @param init_wfc initial wavefunction type [in]
     * @param nbands number of bands [in]
     */
    void cal_nwfc(std::ofstream& log, UnitCell& ucell,Atom* atoms, const int nspin, const int nlocal, const int npol,
               const std::string& basis_type, const std::string& esolver_type, const std::string& init_wfc, const int nbands);

    /**
     * @brief Calculate meshx.
     *
     * Demand: atoms[].msh
     *
     * @param meshx output mesh size [out]
     * @param atoms atom pointer [in]
     * @param ntype number of atom types [in]
     */
    void cal_meshx(int& meshx,const Atom* atoms, const int ntype);

    /**
     * @brief Calculate natomwfc.
     *
     * Demand: atoms[].nchi, atoms[].lchi, atoms[].oc, atoms[].na
     *
     * @param log output file stream [in]
     * @param natomwfc output number of atomic wavefunctions [out]
     * @param ntype number of atom types [in]
     * @param atoms atom pointer [in]
     * @param nspin number of spin components [in]
     */
    void cal_natomwfc(std::ofstream& log,int& natomwfc,const int ntype,const Atom* atoms,const int nspin);

    /**
     * @brief Check consistency between two atom labels from STRU and pseudo or
     * orb file.
     *
     * @param label1 atom label from STRU [in]
     * @param label2 atom label from pseudo or orbital file [in]
     */
    void compare_atom_labels(const std::string& label1, const std::string& label2);

}

#endif // READ_PP_UCELL_H