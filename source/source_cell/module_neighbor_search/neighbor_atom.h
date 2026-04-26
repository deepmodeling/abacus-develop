#ifndef NEIGHBOR_ATOM_H
#define NEIGHBOR_ATOM_H

#include <vector>

class NeighborAtom
{
public:
    double position_x;
    double position_y;
    double position_z;
    int atom_type;
    int atom_index;
    int atom_id;
    //bool isghost;
    bool is_inside;

    NeighborAtom(double x, double y, double z, int type, int index, int id)
        : position_x(x), position_y(y), position_z(z),
          atom_type(type), atom_index(index), atom_id(id) {}
};

class InputAtoms
{
public:
    std::vector<NeighborAtom> InputAtom;
    double x_low, x_high, y_low, y_high, z_low, z_high;
    int n_atoms;

    InputAtoms()
        : x_low(0), x_high(0), y_low(0), y_high(0), z_low(0), z_high(0), n_atoms(0) {}
};

#endif // NEIGHBOR_ATOM_H