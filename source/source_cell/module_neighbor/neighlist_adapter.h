#ifndef NEIGHLIST_ADAPTER_H
#define NEIGHLIST_ADAPTER_H

#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/module_neighlist/neighbor_search.h"

void convert_neighbor_search_to_adjs(const UnitCell& ucell,
                                     const NeighborSearch& neighbor_search,
                                     int ntype,
                                     int nnumber,
                                     AdjacentAtomInfo& adjs);

#endif
