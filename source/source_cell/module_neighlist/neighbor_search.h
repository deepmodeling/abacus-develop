#ifndef NEIGHBOR_SEARCH_H
#define NEIGHBOR_SEARCH_H


#include "source_cell/module_neighlist/neighbor_atom.h"
#include "source_cell/module_neighlist/bin_manager.h"
#include "source_cell/module_neighlist/neighbor_list.h"
#include "source_cell/module_neighlist/unitcell_interface.h"

class NeighborSearch
{
public:
    NeighborSearch()=default;
    ~NeighborSearch()=default;

    void init(const IAtomProvider& ucell, double sr, int mpi_rank);


    void build_neighbors();

    InputAtoms ucell_to_input_atoms(const IAtomProvider& ucell);

    void Check_Expand_Condition(const IAtomProvider& ucell);

    void setMemberVariables(const IAtomProvider& ucell);
    NeighborList& get_neighbor_list() { return neighbor_list; }

    double distance(
        double position_x,
        double position_y,
        double position_z,
        double x_low,
        double y_low,
        double z_low);

    void decompose(int mpi_size, int& nx, int& ny, int& nz);

    double search_radius;


    int x, y, z;
    double wide_x, wide_y, wide_z;

    int glayerX, glayerY, glayerZ;
    int glayerX_minus, glayerY_minus, glayerZ_minus;

    std::vector<NeighborAtom> all_atoms;
    std::vector<NeighborAtom> inside_atoms;
    std::vector<NeighborAtom> ghost_atoms;

    NeighborList neighbor_list;
    BinManager bin_manager;
};

#endif