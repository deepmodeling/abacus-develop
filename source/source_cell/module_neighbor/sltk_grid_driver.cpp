#include "sltk_grid_driver.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#ifdef _OPENMP
#include <omp.h>
#endif

Grid_Driver::Grid_Driver(
	const int &test_d_in, 
	const int &test_grid_in)
:test_deconstructor(test_d_in),
Grid(test_grid_in)
{
	test_deconstructor	= test_d_in;
}

Grid_Driver::~Grid_Driver()
{
}

void Grid_Driver::Find_atom(const UnitCell& ucell,
                            const int ntype,
                            const int nnumber,
                            AdjacentAtomInfo* adjs) const
{
    ModuleBase::timer::start("Grid_Driver", "Find_atom");
    // The public interface now uses the Phase 2.1 AtomPack path by default.
    // Find_atom_from_legacy() is kept only as a regression and fallback route.
    this->Find_atom_from_atom_pack(ucell, ntype, nnumber, adjs);
    ModuleBase::timer::end("Grid_Driver", "Find_atom");
    return;
}

void Grid_Driver::Find_atom_from_legacy(const UnitCell& ucell,
                                        const int ntype,
                                        const int nnumber,
                                        AdjacentAtomInfo* adjs) const
{
    //	std::cout << "lenght in Find atom = " << atomlink[offset].fatom.getAdjacentSet()->getLength() << std::endl;

    // store result in member adj_info when parameter adjs is NULL
    AdjacentAtomInfo* local_adjs = adjs == nullptr ? &this->adj_info : adjs;
    local_adjs->clear();
    if (ntype < 0 || ntype >= static_cast<int>(all_adj_info.size())
        || nnumber < 0 || nnumber >= static_cast<int>(all_adj_info[ntype].size()))
    {
        throw std::runtime_error("Legacy Grid neighbor path is not built for this atom.");
    }

    const std::vector<FAtom*>& all_atom = all_adj_info[ntype][nnumber];

    for (const FAtom* atom: all_atom)
    {
        local_adjs->ntype.push_back(atom->type);
        local_adjs->natom.push_back(atom->natom);
        local_adjs->box.push_back(ModuleBase::Vector3<int>(atom->cell_x, atom->cell_y, atom->cell_z));
        local_adjs->adjacent_tau.push_back(ModuleBase::Vector3<double>(atom->x, atom->y, atom->z));
        local_adjs->adj_num++;
    }
    // 20241204 zhanghaochong
    // for some unknown reason, the last neighbour atom must be it self
    // is self must in last, the order cannot be changed.
    // if self not in last, test 701_LJ_MD_Anderson will assert
	local_adjs->ntype.push_back(ntype);
	local_adjs->natom.push_back(nnumber);
	local_adjs->box.push_back(ModuleBase::Vector3<int>(0, 0, 0));
	local_adjs->adjacent_tau.push_back(ModuleBase::Vector3<double>(ucell.atoms[ntype].tau[nnumber].x, ucell.atoms[ntype].tau[nnumber].y, ucell.atoms[ntype].tau[nnumber].z));
    return;
}

void Grid_Driver::Find_atom_from_atom_pack(const UnitCell& ucell,
                                           const int ntype,
                                           const int nnumber,
                                           AdjacentAtomInfo* adjs) const
{
    AdjacentAtomInfo* local_adjs = adjs == nullptr ? &this->adj_info : adjs;
    local_adjs->clear();
    if (ntype < 0 || ntype >= static_cast<int>(neighbor_pair_indices.size())
        || nnumber < 0 || nnumber >= static_cast<int>(neighbor_pair_indices[ntype].size()))
    {
        throw std::runtime_error("AtomPack neighbor path is not built for this atom.");
    }

    const std::vector<int>& pair_indices = neighbor_pair_indices[ntype][nnumber];
    for (const int pair_index: pair_indices)
    {
        const ModuleNeighbor::NeighborPair& pair = neighbor_pairs[pair_index];
        const int atom_index = pair.neighbor_index;
        if (atom_index < 0 || atom_index >= atom_pack.size())
        {
            throw std::runtime_error("AtomPack neighbor index is out of range.");
        }

        local_adjs->ntype.push_back(pair.neighbor_type);
        local_adjs->natom.push_back(pair.neighbor_natom);
        local_adjs->box.push_back(ModuleBase::Vector3<int>(pair.cell_x, pair.cell_y, pair.cell_z));
        local_adjs->adjacent_tau.push_back(
            ModuleBase::Vector3<double>(atom_pack.x[atom_index], atom_pack.y[atom_index], atom_pack.z[atom_index]));
        local_adjs->adj_num++;
    }

    // Keep the ABACUS compatibility rule from the legacy path: the center atom
    // itself is appended after all real neighbors, and adj_num counts only the
    // real neighbor entries.
    local_adjs->ntype.push_back(ntype);
    local_adjs->natom.push_back(nnumber);
    local_adjs->box.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    local_adjs->adjacent_tau.push_back(ModuleBase::Vector3<double>(ucell.atoms[ntype].tau[nnumber].x,
                                                                   ucell.atoms[ntype].tau[nnumber].y,
                                                                   ucell.atoms[ntype].tau[nnumber].z));
}

void Grid_Driver::Find_atom(const UnitCell& ucell,
                   const ModuleBase::Vector3<double>& cartesian_posi,
                   const int& ntype,
                   const int& nnumber,
                   AdjacentAtomInfo* adjs) const
{
    this->Find_atom(ucell, ntype, nnumber, adjs);
}

// filter_adjs delete not adjacent atoms in adjs
void filter_adjs(const std::vector<bool>& is_adj, AdjacentAtomInfo& adjs)
{
    const int size = adjs.adj_num + 1;
    for (int i = size - 1; i >= 0; --i)
    {
        if (!is_adj[i])
        {
            adjs.adj_num--;
            adjs.ntype.erase(adjs.ntype.begin() + i);
            adjs.natom.erase(adjs.natom.begin() + i);
            adjs.adjacent_tau.erase(adjs.adjacent_tau.begin() + i); // info of adjacent_tau is not used in future
            adjs.box.erase(adjs.box.begin() + i);
        }
    }
}
