#ifndef CMD_NEIGHBOR_H
#define CMD_NEIGHBOR_H

#include "../module_base/vector3.h"
#include "../module_cell/unitcell_pseudo.h"

class CMD_neighbor
{
public:

    CMD_neighbor();
    ~CMD_neighbor();

    ModuleBase::Vector3<double> cell_periodic(const ModuleBase::Vector3<double> a, const ModuleBase::Vector3<double> b);
    void neighbor(const UnitCell_pseudo &ucell_c);
    void comm_list(const int num, int *nlist_in, int **list_in, int *nlist_out, int **list_out);

    int **list;    // record the index of adjent atoms of every atom
    int *nlist;    // record the adjent num of every atom
    int dim;
};
#endif