#ifndef TOQO_H
#define TOQO_H

#include <iostream>
#include <string>
#include "module_cell/unitcell.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_base/atom_in.h"

class toQO
{
    public:
        toQO(std::string qo_basis);
        ~toQO();
        // ----
        // main functions, implemented in to_qo.cpp
        // ----
        /// @brief initialize the QO class
        /// @param p_ucell interface (raw pointer) to the unitcell
        /// @param nkpts number of kpoints
        void initialize(UnitCell *p_ucell, 
                        int nkpts);
        /// @brief calculate the overlap between atomic orbitals and numerical atomic orbitals
        void cal_ovlp_ao_nao();
        /// @brief write the overlap to file
        void write_ovlp_ao_nao();

        // ----
        // tool functions, implemented in to_qo_tools.cpp
        // ----
        /// @brief calculate vectors connecting all atom pairs that needed to calculate their overlap
        void cal_two_center_vectors();

        /// @brief for one certain atom, with given cutoff radius of orbital, check in which supercell the tail of the orbital will vanish
        /// @details given coordiante of one atom, translate it to the supercell again-and-again, until the tail of the orbital vanishes. Finally this function will yield the maximum supercell size
        void scan_supercell_for_atom(int it, int ia, double rcut);

        /// @brief once atom-by-atom supercell scan is done, get maximum supercell size and return
        void scan_supercell();

    private:
        /// @brief interface to the unitcell
        UnitCell *p_ucell_ = nullptr;
        /// @brief number of kpoints
        int nkpts_ = 0;
        /// @brief current atomic orbital basis for generating QO
        std::string qo_basis_ = "hydrogen";

        std::unique_ptr<RadialCollection> nao_;
        std::unique_ptr<RadialCollection> ao_;

        atom_in atom_database_;
};
#endif // TOQO_H