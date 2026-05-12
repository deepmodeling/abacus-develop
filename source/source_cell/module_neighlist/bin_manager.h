#ifndef BIN_MANAGER_H
#define BIN_MANAGER_H

#include <vector>
#include "source_cell/module_neighlist/neighbor_atom.h"
#include "source_cell/module_neighlist/neighbor_list.h"



class Bin
{
public:
    Bin()=default;
    ~Bin()=default;

    int id_x;
    int id_y;
    int id_z;

    std::vector<NeighborAtom> atoms;
};


class BinManager
{
public:

    void init_bins(
        double sr,
        const std::vector<NeighborAtom>& inside_atoms,
        const std::vector<NeighborAtom>& ghost_atoms
    );


    void do_binning(
        const std::vector<NeighborAtom>& inside_atoms,
        const std::vector<NeighborAtom>& ghost_atoms
    );


    void build_atom_neighbors(
        NeighborList& neighbor_list,
        std::vector<NeighborAtom>& atoms
    );


    void clear();


    double sradius;

    double x_min, y_min, z_min;
    double x_max, y_max, z_max;

    double bin_sizex;
    double bin_sizey;
    double bin_sizez;

    int nbinx;
    int nbiny;
    int nbinz;

    std::vector<Bin> bins;
};

#endif // BIN_MANAGER_H