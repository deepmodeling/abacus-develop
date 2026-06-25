#include "sltk_grid_driver.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>

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
    if (mpi_mode() && !is_local_center(ntype, nnumber))
    {
        throw std::runtime_error("MPI Grid neighbor query requires a local center atom.");
    }
    ModuleBase::timer::start("Grid_Driver", "Find_atom");
    // The integer-indexed AtomPack list is the default query path.
    // Find_atom_from_legacy() remains available for regression comparisons.
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
    if (mpi_mode() && !is_local_center(ntype, nnumber))
    {
        throw std::runtime_error("MPI Grid neighbor query requires a local center atom.");
    }
    if (paged_neighbor_list.type_offset.empty())
    {
        throw std::runtime_error("AtomPack neighbor path is not built for this atom.");
    }

    const ModuleBase::Vector3<double>& center_tau = ucell.atoms[ntype].tau[nnumber];
    ModuleBase::Matrix3 inverse_lattice;
    ModuleBase::Vector3<double> center_winding;
    bool center_winding_initialized = false;
    if (pbc)
    {
        inverse_lattice = ucell.latvec.Inverse();
    }

    paged_neighbor_list.for_each_pair_index(
        ntype,
        nnumber,
        [&](const int pair_index)
        {
            if (pair_index < 0 || pair_index >= static_cast<int>(neighbor_pairs.size()))
            {
                throw std::runtime_error("Paged neighbor pair index is out of range.");
            }
            const ModuleNeighbor::NeighborPair& pair = neighbor_pairs[pair_index];
            const int atom_index = pair.neighbor_index;
            if (pair.center_index < 0 || pair.center_index >= atom_pack.size()
                || atom_index < 0 || atom_index >= atom_pack.size())
            {
                throw std::runtime_error("AtomPack pair index is out of range.");
            }

            if (pair.neighbor_type < 0 || pair.neighbor_type >= ucell.ntype
                || pair.neighbor_natom < 0
                || pair.neighbor_natom >= ucell.atoms[pair.neighbor_type].na)
            {
                throw std::runtime_error("Neighbor pair refers to an invalid UnitCell atom.");
            }

            // Pair identities and cell shifts remain valid while the Verlet
            // list is reused, but AtomPack coordinates are rebuild-time
            // snapshots. Reconstruct the current image coordinate so moving
            // atoms never expose stale positions to downstream consumers.
            const ModuleBase::Vector3<double>& base_tau
                = ucell.atoms[pair.neighbor_type].tau[pair.neighbor_natom];
            int cell_x = pair.cell_x;
            int cell_y = pair.cell_y;
            int cell_z = pair.cell_z;
            if (pbc)
            {
                // MD/relax may wrap a physical atom back into the primary cell
                // without rebuilding the Verlet list. Preserve the pair's
                // continuous image by compensating both atoms' integer winding
                // since the last rebuild.
                if (!center_winding_initialized)
                {
                    const ModuleBase::Vector3<double> old_center(atom_pack.x[pair.center_index],
                                                                  atom_pack.y[pair.center_index],
                                                                  atom_pack.z[pair.center_index]);
                    center_winding = old_center * inverse_lattice - center_tau * inverse_lattice;
                    center_winding_initialized = true;
                }
                const ModuleBase::Vector3<double> old_neighbor_image(atom_pack.x[atom_index],
                                                                     atom_pack.y[atom_index],
                                                                     atom_pack.z[atom_index]);
                const ModuleBase::Vector3<double> old_neighbor
                    = old_neighbor_image - static_cast<double>(pair.cell_x) * ucell.a1
                      - static_cast<double>(pair.cell_y) * ucell.a2
                      - static_cast<double>(pair.cell_z) * ucell.a3;
                const ModuleBase::Vector3<double> neighbor_winding
                    = old_neighbor * inverse_lattice - base_tau * inverse_lattice;
                cell_x += static_cast<int>(std::round(neighbor_winding.x))
                          - static_cast<int>(std::round(center_winding.x));
                cell_y += static_cast<int>(std::round(neighbor_winding.y))
                          - static_cast<int>(std::round(center_winding.y));
                cell_z += static_cast<int>(std::round(neighbor_winding.z))
                          - static_cast<int>(std::round(center_winding.z));
            }
            const ModuleBase::Vector3<double> neighbor_tau
                = base_tau + static_cast<double>(cell_x) * ucell.a1
                  + static_cast<double>(cell_y) * ucell.a2
                  + static_cast<double>(cell_z) * ucell.a3;
            const ModuleBase::Vector3<double> delta = neighbor_tau - center_tau;
            const double distance2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
            if (distance2 == 0.0 || distance2 > query_radius2)
            {
                return;
            }

            local_adjs->ntype.push_back(pair.neighbor_type);
            local_adjs->natom.push_back(pair.neighbor_natom);
            local_adjs->box.push_back(ModuleBase::Vector3<int>(cell_x, cell_y, cell_z));
            local_adjs->adjacent_tau.push_back(neighbor_tau);
            local_adjs->adj_num++;
        });

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
