#include "source_cell/module_neighbor/neighlist_adapter.h"

#include <cassert>

void convert_neighbor_search_to_adjs(const UnitCell& ucell,
                                     const NeighborSearch& neighbor_search,
                                     int ntype,
                                     int nnumber,
                                     AdjacentAtomInfo& adjs)
{
    adjs.clear();

    int center_index = -1;
    for (int i = 0; i < static_cast<int>(neighbor_search.inside_atoms.size()); ++i)
    {
        const NeighborAtom& atom = neighbor_search.inside_atoms[i];
        if (atom.atom_type == ntype && atom.atom_index == nnumber && atom.cell_x == 0 && atom.cell_y == 0
            && atom.cell_z == 0)
        {
            center_index = i;
            break;
        }
    }

    assert(center_index >= 0);
    assert(center_index < static_cast<int>(neighbor_search.neighbor_list.numneigh.size()));

    if (center_index >= 0 && center_index < static_cast<int>(neighbor_search.neighbor_list.numneigh.size()))
    {
        const int neighbor_count = neighbor_search.neighbor_list.numneigh[center_index];
        const int* neighbor_ids = neighbor_search.neighbor_list.firstneigh[center_index];
        assert(neighbor_ids != nullptr || neighbor_count == 0);

        for (int i = 0; neighbor_ids != nullptr && i < neighbor_count; ++i)
        {
            const int atom_id = neighbor_ids[i];
            assert(atom_id >= 0);
            assert(atom_id < static_cast<int>(neighbor_search.all_atoms.size()));

            const NeighborAtom& atom = neighbor_search.all_atoms[atom_id];
            adjs.ntype.push_back(atom.atom_type);
            adjs.natom.push_back(atom.atom_index);
            adjs.box.push_back(ModuleBase::Vector3<int>(atom.cell_x, atom.cell_y, atom.cell_z));
            adjs.adjacent_tau.push_back(
                ModuleBase::Vector3<double>(atom.position_x, atom.position_y, atom.position_z));
            adjs.adj_num++;
        }
    }

    adjs.ntype.push_back(ntype);
    adjs.natom.push_back(nnumber);
    adjs.box.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    adjs.adjacent_tau.push_back(ModuleBase::Vector3<double>(ucell.atoms[ntype].tau[nnumber].x,
                                                            ucell.atoms[ntype].tau[nnumber].y,
                                                            ucell.atoms[ntype].tau[nnumber].z));
}
