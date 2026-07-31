/**
 * @file atom_spec.h
 * @brief Atom class for storing atom information.
 */
#ifndef ATOM_H
#define ATOM_H

#include "atom_pseudo.h"

/**
 * @brief Atom class for storing atom information.
 */
class Atom
{
  public:
    /**
     * @brief Constructor.
     */
    Atom();

    /**
     * @brief Destructor.
     */
    ~Atom();

    Atom_pseudo ncpp;                             ///< pseudopotential for this atom type
    double mass = 0.0;                            ///< the mass of atom
    std::vector<ModuleBase::Vector3<int>> mbl;    ///< whether the atoms can move or not
    bool flag_empty_element = false;              ///< whether is the empty element for bsse (Peize Lin add 2021.04.07)

    std::vector<int> iw2m;                       ///< use iw to find m
    std::vector<int> iw2n;                       ///< use iw to find n
    std::vector<int> iw2l;                       ///< use iw to find L
    std::vector<int> iw2_ylm;                    ///< use iw to find ylm index
    std::vector<bool> iw2_new;                   ///< use iw to find new flag
    int nw = 0;                                  ///< number of local orbitals (l,n,m) of this type

    /**
     * @brief Set index arrays.
     */
    void set_index();

    int type = 0;                                ///< Index of atom type
    int na = 0;                                  ///< Number of atoms in this type

    int nwl = 0;                                 ///< max L(Angular momentum) (for local basis)
    double Rcut = 0.0;                           ///< cut-off radius (pengfei Li 16-2-29)
    std::vector<int> l_nchi;                     ///< number of chi for each L
    int stapos_wf = 0;                           ///< start position of wave functions

    std::string label;                            ///< atomic symbol
    std::vector<ModuleBase::Vector3<double>> tau; ///< Cartesian coordinates of each atom in this type
    std::vector<ModuleBase::Vector3<double>> dis; ///< direct displacements of each atom in this type in current step (liuyu modify 2023-03-22)
    std::vector<ModuleBase::Vector3<double>> taud; ///< Direct coordinates of each atom in this type
    std::vector<ModuleBase::Vector3<int>> boundary_shift; ///< record for periodic boundary adjustment
    std::vector<ModuleBase::Vector3<double>> vel; ///< velocities of each atom in this type
    std::vector<ModuleBase::Vector3<double>> force; ///< force acting on each atom in this type
    std::vector<ModuleBase::Vector3<double>> lambda; ///< Lagrange multiplier for each atom in this type, used in deltaspin
    std::vector<ModuleBase::Vector3<int>> constrain; ///< constrain for each atom in this type, used in deltaspin
    std::string label_orb;                          ///< atomic element symbol in the orbital file of lcao

    std::vector<double> mag;                      ///< magnetic moment
    std::vector<double> angle1;                   ///< spin angle, added by zhengdy-soc
    std::vector<double> angle2;                   ///< spin angle, added by zhengdy-soc
    std::vector<ModuleBase::Vector3<double>> m_loc_; ///< local magnetic moment

    /// @brief Coulomb potential v(r) = z/r
    /// It is a local potential, and has no non-local potential parts.
    bool coulomb_potential = false;

    /**
     * @brief Print atom information.
     *
     * @param ofs output file stream
     */
    void print_Atom(std::ofstream& ofs);

    /**
     * @brief Update force.
     *
     * @param fcs force constant matrix
     */
    void update_force(ModuleBase::matrix& fcs);

#ifdef __MPI
    /**
     * @brief Broadcast atom data.
     */
    void bcast_atom();

    /**
     * @brief Broadcast atom data (second version).
     */
    void bcast_atom2();
#endif
};

#endif // Atomspec
