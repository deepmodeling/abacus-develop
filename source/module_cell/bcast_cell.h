#ifndef BCAST_CELL_H
#define BCAST_CELL_H

#include "module_cell/unitcell.h"
namespace unitcell
{
    /**
     * @brief broadcast the tau array of the atoms
     * 
     * @param atoms: the atoms to be broadcasted [in/out]
     * @param ntype: the number of types of the atoms [in]
     */
    void bcast_atoms_tau(Atom* atoms,
                     const int ntype);
                     
    /**
     * @brief broadcast the pseduo of the atoms
     * 
     * @param atoms: the atoms to be broadcasted [in/out]
     * @param ntype:  the number of types of the atoms [in]
     */
    void bcast_atoms_pseudo(Atom* atoms,
                            const int ntype);
    /**
     * @brief broadcast the lattice
     * 
     * @param lat: the lattice to be broadcasted [in/out]
     */
    void bcast_Lattice(Lattice& lat);

    /**
     * @brief broadcast the magnetism
     * 
     * @param magnet: the magnetism to be broadcasted [in/out]
     * @param nytpe: the number of types of the atoms [in]
     */
    void bcast_magnetism(Magnetism& magnet,
                         const int ntype);
    
    /**
     * @brief broadcast the unitcell
     * 
     * @param ucell: the unitcell to be broadcasted [in/out]
     */
    void bcast_unitcell(UnitCell& ucell);


}

#endif // BCAST_CELL_H